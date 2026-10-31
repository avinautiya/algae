import copy
import itertools
import json
import math
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ice_incision_forecast.core import ThermalMap,ConductancePrior,region_bounds
from ice_incision_forecast.forecast import forecast_case
from ice_incision_forecast.__main__ import demo_case
from ice_incision_forecast.validation import validate_randomized

ROOT=Path(__file__).resolve().parents[2]


def rk4(model,conductance,steps):
    p=model.inlet_power_W; q=model.inlet_Q_m3_s; dg=conductance/steps
    def rhs(p,q):
        power=p/(model.rho_water*model.cp*q)
        return -power,power/(model.rho_water*model.latent_heat)
    for _ in range(steps):
        k1=rhs(p,q); k2=rhs(p+dg*k1[0]/2,q+dg*k1[1]/2)
        k3=rhs(p+dg*k2[0]/2,q+dg*k2[1]/2); k4=rhs(p+dg*k3[0],q+dg*k3[1])
        p+=dg*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        q+=dg*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
    return p,q


class ThermalMapTests(unittest.TestCase):
    def test_mass_and_energy_invariant(self):
        model=ThermalMap(0.02,0.3)
        for g in (0,1000,10000,100000,1e7):
            r=model.state(g)
            self.assertAlmostEqual(r["power_W"]+r["removed_W"],model.inlet_power_W,places=9)
            self.assertAlmostEqual((r["Q_m3_s"]-model.inlet_Q_m3_s)*model.rho_water*model.latent_heat,r["removed_W"],delta=1e-8)

    def test_independent_rk4_fourth_order(self):
        model=ThermalMap(0.02,0.3); g=100000; exact=model.state(g)
        errors=[]
        for n in (16,32,64):
            p,q=rk4(model,g,n)
            errors.append(abs(p-exact["power_W"])/model.inlet_power_W)
            self.assertAlmostEqual(q,exact["Q_m3_s"],delta=2e-10)
        self.assertGreater(errors[0]/errors[1],12)
        self.assertGreater(errors[1]/errors[2],12)
        self.assertLess(errors[-1],1e-8)

    def test_attenuation_bracket_and_convexity(self):
        model=ThermalMap(0.01,2)
        for g in (1,1000,10000):
            s=model.attenuation(g); target=g/model.capacity_W_K
            self.assertLessEqual(target/(1+model.ste)-1e-14,s)
            self.assertLessEqual(s,target+1e-14)
        values=[model.state(g)["power_W"] for g in (1000,2000,3000)]
        self.assertGreater(values[0]-2*values[1]+values[2],0)

    def test_temperature_calibration_roundtrip(self):
        model=ThermalMap(0.1,0.2)
        for g in (0,1000,20000,100000):
            t=model.state(g)["temperature_excess_K"]
            self.assertAlmostEqual(model.total_conductance_from_temperature(t),g,delta=1e-8)

    def test_zero_heat_and_mass_addition_limit(self):
        model=ThermalMap(0.1,0)
        self.assertEqual(model.state(1e6)["Q_m3_s"],0.1)
        self.assertEqual(model.region_power(0,1e6),0)
        model=ThermalMap(0.1,1e-8); g=1e4
        self.assertAlmostEqual(model.state(g)["power_W"]/model.inlet_power_W,math.exp(-g/model.capacity_W_K),places=10)

    def test_invalid_inputs_and_zero_outlet_calibration(self):
        with self.assertRaises(ValueError):
            ThermalMap(0,0.2)
        with self.assertRaises(ValueError):
            ThermalMap(.1,.2).total_conductance_from_temperature(0)
        with self.assertRaises(ValueError):
            ThermalMap(.1,.2).state(-1)


class SharpBoundTests(unittest.TestCase):
    def test_continuous_polytope_reference(self):
        result=validate_randomized(cases=30,seed=991)
        self.assertLess(result["max_power_scaled_error"],1e-10)

    def exhaustive(self,prior,units=10000):
        n=len(prior.lower_W_K)
        candidates=[list(x) for x in itertools.product(*[range(round(l/units),round(u/units)+1) for l,u in zip(prior.lower_W_K,prior.upper_W_K)])]
        candidates=[[units*x for x in p] for p in candidates if prior.feasible([units*x for x in p])]
        self.assertTrue(candidates)
        model=ThermalMap(.02,.3)
        for a in range(n):
            for b in range(a+1,n+1):
                r=region_bounds(prior,a,b)
                values=[model.region_power(sum(x[:a]),sum(x[a:b])) for x in candidates]
                for label,reference in (("minimum",min(values)),("maximum",max(values))):
                    w=r[label]
                    self.assertTrue(prior.feasible(w["profile_W_K"]))
                    predicted=model.region_power(w["prefix_W_K"],w["region_W_K"])
                    self.assertAlmostEqual(predicted,reference,delta=1e-8)

    def test_exact_total_exhaustive_vertices(self):
        self.exhaustive(ConductancePrior((0,)*4,(30000,)*4,50000,50000))

    def test_interval_total_exhaustive_vertices(self):
        self.exhaustive(ConductancePrior((10000,0,10000,0),(30000,20000,30000,30000),50000,70000))

    def test_exact_internal_sensors_exhaustive_vertices(self):
        self.exhaustive(ConductancePrior((0,)*6,(30000,)*6,90000,90000,((2,30000),(4,60000))))

    def test_internal_sensors_and_final_total_band(self):
        self.exhaustive(ConductancePrior((0,)*4,(30000,)*4,40000,70000,((2,30000),)))

    def test_one_fixed_witness_attains_time_integrated_bound(self):
        prior=ConductancePrior((0,)*4,(30000,)*4,50000,50000)
        candidates=[[10000*x for x in p] for p in itertools.product(range(4),repeat=4) if sum(p)==5]
        scenarios=[(ThermalMap(.01,.1),300,0.8),(ThermalMap(.1,.4),900,1.4),(ThermalMap(.03,0),500,1.0)]
        r=region_bounds(prior,1,3)
        def total(profile):
            return sum(m.region_power(mult*sum(profile[:1]),mult*sum(profile[1:3]))*dt for m,dt,mult in scenarios)
        values=[total(p) for p in candidates]
        self.assertAlmostEqual(total(r["minimum"]["profile_W_K"]),min(values),delta=1e-5)
        self.assertAlmostEqual(total(r["maximum"]["profile_W_K"]),max(values),delta=1e-5)

    def test_outlet_equivalence_does_not_fix_regional_melt(self):
        prior=ConductancePrior((200,)*20,(1800,)*20,20000,20000)
        r=region_bounds(prior,8,12)
        a,b=r["minimum"],r["maximum"]
        for q,t in ((.02,.1),(.1,.3),(.04,.2)):
            model=ThermalMap(q,t)
            self.assertEqual(model.state(a["total_W_K"]),model.state(b["total_W_K"]))
            self.assertGreater(model.region_power(b["prefix_W_K"],b["region_W_K"]),
                               5*model.region_power(a["prefix_W_K"],a["region_W_K"]))

    def test_bracketing_exact_sensors_identify_regional_heat(self):
        prior=ConductancePrior((200,)*20,(1800,)*20,20000,20000,((8,8000),(12,12000)))
        r=region_bounds(prior,8,12)
        self.assertAlmostEqual(r["minimum"]["prefix_W_K"],r["maximum"]["prefix_W_K"])
        self.assertAlmostEqual(r["minimum"]["region_W_K"],r["maximum"]["region_W_K"])

    def test_randomized_feasible_transfers_stay_within_bounds(self):
        rng=random.Random(771)
        prior=ConductancePrior((200,)*10,(1800,)*10,10000,10000)
        model=ThermalMap(.02,.2); r=region_bounds(prior,3,6)
        lo=model.region_power(r["minimum"]["prefix_W_K"],r["minimum"]["region_W_K"])
        hi=model.region_power(r["maximum"]["prefix_W_K"],r["maximum"]["region_W_K"])
        p=[1000.0]*10
        for _ in range(1000):
            i,j=rng.sample(range(10),2)
            delta=rng.uniform(max(200-p[i],p[j]-1800),min(1800-p[i],p[j]-200))
            p[i]+=delta; p[j]-=delta
            value=model.region_power(sum(p[:3]),sum(p[3:6]))
            self.assertGreaterEqual(value,lo-1e-8); self.assertLessEqual(value,hi+1e-8)

    def test_infeasible_priors_and_invalid_regions(self):
        with self.assertRaises(ValueError):
            ConductancePrior((1,1),(2,2),5,6)
        with self.assertRaises(ValueError):
            ConductancePrior((1,1),(2,2),3,3,((1,.5),))
        p=ConductancePrior((1,1),(2,2),3,3)
        with self.assertRaises(ValueError):
            region_bounds(p,1,1)


class ForecastTests(unittest.TestCase):
    def test_complete_forward_envelope_and_ledgers(self):
        result=forecast_case(demo_case())
        self.assertLess(result["summary"]["max_abs_water_residual_m3"],1e-10)
        self.assertLess(result["summary"]["max_abs_energy_residual_J"],1e-6)
        self.assertFalse(result["summary"]["scope_flags"])
        for r in result["bed"]:
            self.assertLessEqual(r["lowering_min_m"]-1e-12,r["lowering_nominal_m"])
            self.assertLessEqual(r["lowering_nominal_m"],r["lowering_max_m"]+1e-12)

    def test_fixed_geometry_scope_violation_is_flagged(self):
        case=demo_case(); case["scope_checks"]["max_fractional_depth_change"]=1e-6
        self.assertIn("FROZEN_GEOMETRY_CHANGE_THRESHOLD_EXCEEDED",forecast_case(case)["summary"]["scope_flags"])

    def test_minimax_forecast_and_attainable_error_radius(self):
        result=forecast_case(demo_case())
        for r in result["regions"]:
            lower,upper=r["lowering_min_m"],r["lowering_max_m"]
            center,radius=r["lowering_minimax_m"],r["minimax_worst_case_error_m"]
            self.assertAlmostEqual(max(abs(center-lower),abs(center-upper)),radius,places=14)
            self.assertLessEqual(radius,r["nominal_worst_case_error_m"]+1e-14)
            for offset in (-2*radius,-radius/2,radius/2,2*radius):
                shifted=center+offset
                self.assertGreaterEqual(max(abs(shifted-lower),abs(shifted-upper))+1e-14,radius)

    def test_threshold_occurrence_bounds(self):
        case=demo_case(); r=forecast_case(case)["summary"]["threshold_forecasts"]["declared_control_region"]
        self.assertEqual(r["status"],"OCCURRENCE_UNRESOLVED_WITHIN_HORIZON")
        self.assertGreater(r["earliest_possible_s"],0)
        self.assertIsNone(r["latest_guaranteed_s"])
        case["regions"][0]["critical_lowering_m"]=1e-6
        r=forecast_case(case)["summary"]["threshold_forecasts"]["declared_control_region"]
        self.assertEqual(r["status"],"ALL_PROFILES_CROSS_WITHIN_HORIZON")
        self.assertLessEqual(r["earliest_possible_s"],r["latest_guaranteed_s"])
        case["regions"][0]["critical_lowering_m"]=1.0
        r=forecast_case(case)["summary"]["threshold_forecasts"]["declared_control_region"]
        self.assertEqual(r["status"],"NO_PROFILE_CROSSES_WITHIN_HORIZON")

    def test_future_data_and_invalid_nominal_are_rejected(self):
        case=demo_case(); case["future_outlet_temperature"]=.1
        with self.assertRaises(ValueError):
            forecast_case(case)
        case=demo_case(); case["cells"][0]["conductance_nominal_W_K"]=0
        with self.assertRaises(ValueError):
            forecast_case(case)

    def test_cli_actual_forecast_and_output_preservation(self):
        with tempfile.TemporaryDirectory() as temp:
            case=Path(temp)/"case.json"; out=Path(temp)/"run"
            case.write_text(json.dumps(demo_case()))
            args=[sys.executable,"-m","ice_incision_forecast","predict","--case",str(case),"--output",str(out)]
            result=subprocess.run(args,cwd=ROOT,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            original=(out/"summary.json").read_bytes()
            again=subprocess.run(args,cwd=ROOT,capture_output=True,text=True)
            self.assertNotEqual(again.returncode,0)
            self.assertEqual(original,(out/"summary.json").read_bytes())


if __name__=="__main__":
    unittest.main()
