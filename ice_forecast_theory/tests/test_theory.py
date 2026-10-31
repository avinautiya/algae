import copy
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ice_forecast_theory.forcing import constraint_projector, optimal_forcing, zero_order_hold
from ice_forecast_theory.regime import density, expansion, regime_report

ROOT=Path(__file__).resolve().parents[2]


def scalar_case(rate=-1.0,intervals=2,budget=True):
    dt=2.0/intervals
    return dict(intervals=[dict(duration_s=dt,A_per_s=[[rate]],B_per_s=[[1.0]]) for _ in range(intervals)],
        C=[[1.0]],budget_rows_on_input_amplitudes=[[dt]*intervals] if budget else [],
        state_normalization="Abstract scalar",input_normalization="Abstract scalar",
        output_normalization="Abstract scalar",operator_provenance="Manufactured test only")


class ForcingTests(unittest.TestCase):
    def test_zoh_singular_integrator(self):
        phi,gamma=zero_order_hold([[0]],[[3]],0.4)
        self.assertEqual(phi,[[1]])
        self.assertAlmostEqual(gamma[0][0],1.2,places=14)

    def test_zoh_scalar_against_closed_integral(self):
        phi,gamma=zero_order_hold([[-2]],[[1]],0.7)
        self.assertAlmostEqual(phi[0][0],math.exp(-1.4),places=13)
        self.assertAlmostEqual(gamma[0][0],(1-math.exp(-1.4))/2,places=13)

    def test_unconstrained_integrator_grid_invariant_gain(self):
        for count in (2,4,8):
            r=optimal_forcing(scalar_case(0,count,False))
            self.assertAlmostEqual(r["squared_gain"],4,places=12)
            self.assertAlmostEqual(r["realized_squared_output"],4,places=12)

    def test_equal_budget_integrator_has_zero_terminal_response(self):
        r=optimal_forcing(scalar_case(0))
        self.assertLess(r["squared_gain"],1e-25)
        self.assertLess(abs(r["budget_residuals"][0]),1e-12)

    def test_equal_budget_relaxation_against_closed_formula(self):
        r=optimal_forcing(scalar_case())
        exact=(1-math.exp(-1))**4
        self.assertAlmostEqual(r["squared_gain"],exact,places=13)
        self.assertAlmostEqual(r["realized_squared_output"],exact,places=13)
        amplitudes=[row[0] for row in r["optimal_interval_amplitudes"]]
        self.assertAlmostEqual(sum(amplitudes),0,places=13)
        self.assertAlmostEqual(sum(x*x for x in amplitudes)/2,1,places=13)

    def test_forcing_refinement_below_continuous_time_optimum(self):
        exact=1-math.exp(-4)-(1-math.exp(-2))**2
        values=[optimal_forcing(scalar_case(-1,n))["squared_gain"] for n in (2,4,8)]
        self.assertTrue(values[0]<values[1]<values[2]<exact)

    def test_time_order_noncommuting_operators(self):
        case=scalar_case(budget=False)
        case["intervals"]=[dict(duration_s=1,A_per_s=[[0,1],[0,0]],B_per_s=[[1],[0]]),
                           dict(duration_s=1,A_per_s=[[0,0],[1,0]],B_per_s=[[1],[0]])]
        case["C"]=[[0,1]]
        r=optimal_forcing(case)
        self.assertAlmostEqual(r["squared_gain"],2.5,places=12)
        self.assertAlmostEqual(r["gain_replay_difference"],0,places=12)

    def test_redundant_and_full_rank_budgets(self):
        case=scalar_case()
        case["budget_rows_on_input_amplitudes"]=[[1,1],[2,2],[0,0]]
        r=optimal_forcing(case)
        self.assertEqual(r["budget_rank"],1)
        self.assertAlmostEqual(r["squared_gain"],(1-math.exp(-1))**4,places=12)
        case["budget_rows_on_input_amplitudes"]=[[1,0],[0,1]]
        r=optimal_forcing(case)
        self.assertEqual(r["admissible_coordinate_dimension"],0)
        self.assertEqual(r["squared_gain"],0)

    def test_projector_symmetry_idempotence_and_unit_scale_invariance(self):
        p,rank=constraint_projector([[1e-100,2e-100,3e-100],[2,4,6]],3)
        self.assertEqual(rank,1)
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(p[i][j],p[j][i],places=14)
                self.assertAlmostEqual(sum(p[i][k]*p[k][j] for k in range(3)),p[i][j],places=14)

    def test_invalid_dimensions_norms_and_limits(self):
        with self.assertRaises(ValueError):
            zero_order_hold([[0]],[[1]],0)
        case=scalar_case()
        del case["operator_provenance"]
        with self.assertRaises(ValueError):
            optimal_forcing(case)
        with self.assertRaises(ValueError):
            optimal_forcing(scalar_case(intervals=33))


def regime_case():
    return dict(temperature_C=0.2,depth_m=0.05,velocity_m_s=1,length_m=5,width_m=0.5,
                ice_position="bed",local_patch=dict(depth_m=0.05,speed_m_s=0.005,length_m=0.2))


class RegimeTests(unittest.TestCase):
    def test_published_polynomial_endpoints_and_scale_conversion(self):
        self.assertEqual(density(0),999.842594)
        self.assertAlmostEqual(density(4,"IPTS-68"),999.974958175956,places=9)
        self.assertAlmostEqual(density(3),density(3*1.00024,"IPTS-68"),places=12)

    def test_density_maximum_and_expansion_sign(self):
        self.assertLess(expansion(0.2),0)
        self.assertGreater(expansion(8),0)
        self.assertGreater(density(4),density(3))
        self.assertGreater(density(4),density(5))

    def test_analytic_expansion_against_finite_difference(self):
        t=0.7; dt=0.001
        derivative=(density(t+dt)-density(t-dt))/(2*dt)
        self.assertAlmostEqual(expansion(t),-derivative/density(t),delta=1e-10)

    def test_bed_and_ceiling_reverse_mean_stratification(self):
        case=regime_case(); bed=regime_report(case)["diagnostics"]
        case["ice_position"]="ceiling"; top=regime_report(case)["diagnostics"]
        self.assertLess(bed["buoyancy_N2_per_s2"],0)
        self.assertAlmostEqual(bed["buoyancy_N2_per_s2"],-top["buoyancy_N2_per_s2"],places=12)

    def test_local_slow_flow_changes_inertia_ratio(self):
        r=regime_report(regime_case())["diagnostics"]
        ratio=r["local_patch"]["signed_Ri_N2_h2_over_U2"]/r["signed_Ri_N2_h2_over_U2"]
        self.assertAlmostEqual(ratio,40000,places=8)

    def test_dimensionless_and_power_identities(self):
        case=regime_case(); case["potential_drop_m"]=0.1
        result=regime_report(case); r=result["diagnostics"]; p=result["properties"]
        self.assertAlmostEqual(r["Pe_Uh"],r["Re_Uh"]*r["Pr"],places=8)
        self.assertAlmostEqual(r["We_h"]/r["Bo_h"],r["Fr"]**2,places=12)
        self.assertAlmostEqual(r["potential_to_advected_power_ratio"],p["g_m_s2"]*0.1/(p["cp_J_kg_K"]*0.2),places=14)

    def test_stagnation_and_zero_temperature_do_not_return_infinity(self):
        case=regime_case(); case["velocity_m_s"]=0; case["temperature_C"]=0
        r=regime_report(case)
        json.dumps(r,allow_nan=False)
        self.assertIsNone(r["diagnostics"]["signed_Ri_N2_h2_over_U2"])
        self.assertIsNone(r["diagnostics"]["potential_to_advected_power_ratio"])

    def test_invalid_unsupported_physics_is_rejected(self):
        case=regime_case(); case["salinity"]=35
        with self.assertRaises(ValueError):
            regime_report(case)
        with self.assertRaises(ValueError):
            density(-0.1)

    def test_all_cli_paths_and_immutable_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); demo=root/"demos"
            def run(*args):
                p=subprocess.run([sys.executable,"-m","ice_forecast_theory",*map(str,args)],cwd=ROOT,capture_output=True,text=True)
                self.assertEqual(p.returncode,0,p.stderr)
            run("demo","--output",demo)
            for command,file in (("forcing","forcing_manufactured.json"),("regime","regime_synthetic.json")):
                run(command,"--input",demo/file,"--output",root/(command+".json"))
            p=subprocess.run([sys.executable,"-m","ice_forecast_theory","regime","--input",str(demo/"regime_synthetic.json"),
                              "--output",str(root/"regime.json")],cwd=ROOT,capture_output=True,text=True)
            self.assertNotEqual(p.returncode,0)


if __name__=="__main__":
    unittest.main()
