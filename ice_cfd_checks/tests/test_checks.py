import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ice_cfd_checks.references import DamBreak, StefanReference, rectangular_jump
from ice_cfd_checks.diagnostics import (
    convergence, integrate_faces, norms, score_dam, score_jump, score_rest, score_thermal,
)

ROOT=Path(__file__).resolve().parents[2]


class JumpTests(unittest.TestCase):
    def test_independent_equal_momentum_pair(self):
        # h1=1, h2=2, g=1 requires q^2=3 directly from momentum balance.
        r=rectangular_jump(1,math.sqrt(3),1,gravity=1,rho=1)
        self.assertAlmostEqual(r["depth_down_m"],2,places=14)
        self.assertAlmostEqual(r["head_loss_m"],0.125,places=14)
        self.assertAlmostEqual(r["momentum_up_N"],3.5,places=14)
        self.assertAlmostEqual(r["momentum_down_N"],3.5,places=14)

    def test_jumps_satisfy_momentum_continuity_and_energy_loss(self):
        for h,fr in ((0.01,1.01),(0.2,2),(1,5),(3,12)):
            q=0.7*h*fr*math.sqrt(9.81*h)
            r=rectangular_jump(h,q,0.7)
            self.assertLess(r["Fr_down"],1)
            self.assertAlmostEqual(r["momentum_up_N"]/r["momentum_down_N"],1,places=12)
            self.assertAlmostEqual((r["specific_energy_up_m"]-r["specific_energy_down_m"])/h,
                                   r["head_loss_m"]/h,places=11)

    def test_subcritical_and_nonfinite_rejected(self):
        for h,q,w in ((1,0.1,1),(float("nan"),2,1),(1,2,0)):
            with self.assertRaises(ValueError):
                rectangular_jump(h,q,w)

    def test_scorer_detects_wrong_tailwater(self):
        data=dict(width_m=1,gravity_m_s2=1,rho_kg_m3=1,
            upstream=dict(depth_m=1,discharge_m3_s=math.sqrt(3)),
            downstream=dict(depth_m=2.2,discharge_m3_s=math.sqrt(3)))
        r=score_jump(data)
        self.assertAlmostEqual(r["depth_down_error_relative"],0.1)
        self.assertGreater(abs(r["ideal_momentum_difference_relative"]),0.01)


class DamBreakTests(unittest.TestCase):
    def test_wet_shock_rankine_hugoniot(self):
        r=DamBreak(1,0.1); h,u=r.middle; s=r.speeds[-1]; g=r.gravity
        self.assertAlmostEqual(s*(h-0.1),h*u,places=12)
        self.assertAlmostEqual(s*h*u,h*u*u+g*h*h/2-g*0.1**2/2,places=12)

    def test_rarefaction_invariant(self):
        for r in (DamBreak(1,0.1),DamBreak(1,0)):
            xi=(r.speeds[0]+r.speeds[1])/2
            h,q=r.point(xi*0.2,0.2)
            self.assertAlmostEqual(q/h+2*math.sqrt(r.gravity*h),2*math.sqrt(r.gravity),places=12)

    def test_initial_discontinuity_cell_average(self):
        h,q=DamBreak(1,0.1).cell_average(-0.25,0.75,0)
        self.assertAlmostEqual(h,0.325)
        self.assertEqual(q,0)

    def test_integrated_mass_and_boundary_pressure_impulse(self):
        # Waves remain inside [-10,10]; q=0 at both ends, pressure differs.
        for right in (0,0.1,0.6):
            r=DamBreak(1,right); time=0.3
            h,q=r.cell_average(-10,10,time)
            self.assertAlmostEqual(20*h,10*(1+right),places=11)
            self.assertAlmostEqual(20*q,0.5*r.gravity*(1-right*right)*time,places=11)

    def test_partition_independence_of_cell_integrals(self):
        r=DamBreak(1,0.1)
        a,b,c=-1.2,0.123,1.9
        h,q=r.cell_average(a,c,0.3)
        h0,q0=r.cell_average(a,b,0.3); h1,q1=r.cell_average(b,c,0.3)
        self.assertAlmostEqual(h*(c-a),h0*(b-a)+h1*(c-b),places=12)
        self.assertAlmostEqual(q*(c-a),q0*(b-a)+q1*(c-b),places=12)

    def test_exact_dry_front_and_states(self):
        r=DamBreak(1,0); front=r.speeds[-1]*0.2
        self.assertEqual(r.point(front+0.1,0.2),(0,0))
        self.assertEqual(r.point(-10,0.2),(1,0))
        h,q=r.point(front-1e-6,0.2)
        self.assertGreater(h,0); self.assertGreater(q,0)

    def test_finite_volume_scoring_detects_bias(self):
        r=DamBreak(1,0.1); rows=[]
        for a,b in ((-2,-0.3),(-0.3,0.2),(0.2,2)):
            h,q=r.cell_average(a,b,0.3)
            rows.append(dict(left_m=a,right_m=b,depth_mean_m=h+0.01,hu_mean_m2_s=q))
        result=score_dam(rows,r,0.3)
        self.assertAlmostEqual(result["depth_error_m"]["L1"],0.01)
        self.assertEqual(result["hu_error_m2_s"]["L1"],0)

    def test_invalid_cases_and_overlapping_cells(self):
        with self.assertRaises(ValueError):
            DamBreak(1,2)
        with self.assertRaises(ValueError):
            DamBreak(1,0.1).cell_average(1,0,1)
        rows=[dict(left_m=0,right_m=1,depth_mean_m=1,hu_mean_m2_s=0)]*2
        with self.assertRaises(ValueError):
            score_dam(rows,DamBreak(1,0.1),1)


class StefanTests(unittest.TestCase):
    def test_stefan_condition_and_global_energy(self):
        for ste in (0.001,0.003,0.013,0.2):
            r=StefanReference(ste); time=100
            speed=r.lambda_value*math.sqrt(r.diffusivity/time)
            self.assertAlmostEqual(r.interface_flux(time)/(r.rho_kg_m3*r.latent_heat_J_kg*speed),1,places=12)
            ledger=r.ledger_per_area(time)
            self.assertLess(abs(ledger["residual_J_m2"])/ledger["wall_input_J_m2"],1e-12)

    def test_storage_against_independent_temperature_quadrature(self):
        r=StefanReference(0.013); time=3600; front=r.front(time); n=400
        values=[r.temperature(front*i/n,time) for i in range(n+1)]
        integral=(front/n)/3*(values[0]+values[-1]+4*sum(values[1:-1:2])+2*sum(values[2:-1:2]))
        expected=r.rho_kg_m3*r.cp_J_kg_K*integral
        actual=r.ledger_per_area(time)["liquid_storage_J_m2"]
        self.assertAlmostEqual(actual/expected,1,places=9)

    def test_temperature_satisfies_heat_equation(self):
        r=StefanReference(0.013); time=100; x=r.front(time)*0.4
        dx=r.front(time)*1e-3; dt=0.01
        tder=(r.temperature(x,time+dt)-r.temperature(x,time-dt))/(2*dt)
        xder=(r.temperature(x+dx,time)-2*r.temperature(x,time)+r.temperature(x-dx,time))/(dx*dx)
        self.assertAlmostEqual(tder/(r.diffusivity*xder),1,delta=2e-5)

    def test_boundaries_scaling_and_small_stefan_limit(self):
        r=StefanReference(0.001)
        self.assertEqual(r.front(0),0)
        self.assertAlmostEqual(r.temperature(0,100),r.wall_excess_K)
        self.assertEqual(r.temperature(r.front(100),100),0)
        self.assertAlmostEqual(r.front(400)/r.front(100),2)
        small=StefanReference(1e-10)
        self.assertAlmostEqual(small.lambda_value/math.sqrt(1e-10/2),1,places=9)

    def test_invalid_initial_singular_flux(self):
        with self.assertRaises(ValueError):
            StefanReference(0)
        with self.assertRaises(ValueError):
            StefanReference(0.001).interface_flux(0)


class DiagnosticTests(unittest.TestCase):
    def test_hydrostatic_rest_dry_and_wet(self):
        rows=[dict(bed_m=b,depth_m=max(1-b,0),velocity_m_s=0) for b in (0,0.5,1.2)]
        r=score_rest(rows,1,1)
        self.assertEqual(r["depth_error_m"]["Linf"],0)
        rows[2]["depth_m"]=0.001
        self.assertAlmostEqual(score_rest(rows,1,1)["depth_error_m"]["Linf"],0.001)

    def test_weighted_norms(self):
        r=norms([1,3],[3,1])
        self.assertEqual(r["L1"],1.5)
        self.assertAlmostEqual(r["L2"],math.sqrt(3))

    def test_face_backflow_and_temperature_weighting(self):
        rows=[dict(area_m2=1,liquid_fraction=1,normal_velocity_m_s=2,
                  liquid_volume_flux_m3_s=2,temperature_excess_K=0.2),
              dict(area_m2=1,liquid_fraction=0.5,normal_velocity_m_s=-1,
                  liquid_volume_flux_m3_s=-0.5,temperature_excess_K=0.6)]
        r=integrate_faces(rows)
        self.assertEqual(r["numerical_net_outflow_m3_s"],1.5)
        self.assertAlmostEqual(r["outflow_weighted_temperature_K"],0.2)
        self.assertAlmostEqual(r["inflow_weighted_temperature_K"],0.6)
        self.assertNotEqual(r["wet_area_weighted_temperature_K"],0.2)
        rows[0]["liquid_volume_flux_m3_s"]=0.5
        self.assertEqual(integrate_faces(rows)["numerical_net_outflow_m3_s"],0)

    def test_numerical_flux_is_not_replaced_by_alpha_u(self):
        rows=[dict(area_m2=1,liquid_fraction=0.5,normal_velocity_m_s=2,
              liquid_volume_flux_m3_s=0.8,temperature_excess_K=0.3)]
        r=integrate_faces(rows)
        self.assertEqual(r["numerical_net_outflow_m3_s"],0.8)
        self.assertAlmostEqual(r["reconstructed_minus_numerical_m3_s"],0.2)

    def test_ledger_detects_missing_conductive_flux(self):
        rows=[dict(time_s=t,water_sensible_J=100+2*t,advective_in_J=10*t,
            advective_out_J=5*t,conductive_in_J=t,conductive_out_J=t,
            external_J=t,to_ice_J=4*t) for t in (0,1,2)]
        self.assertEqual(score_thermal(rows)["max_abs_residual_J"],0)
        rows[-1]["conductive_out_J"]=0
        self.assertEqual(score_thermal(rows)["samples"][-1]["residual_J"],-2)

    def test_ledger_restart_origin_and_invalid_time(self):
        rows=[dict(time_s=t,water_sensible_J=100+2*t,advective_in_J=10*t,
            advective_out_J=5*t,conductive_in_J=t,conductive_out_J=t,
            external_J=t,to_ice_J=4*t) for t in (10,11,12)]
        self.assertEqual(score_thermal(rows)["max_abs_residual_J"],0)
        with self.assertRaises(ValueError):
            score_thermal([rows[0],rows[0]])

    def test_convergence_and_oscillation_are_distinguished(self):
        levels=[dict(spacing_m=h,value=2+3*h*h) for h in (0.2,0.1,0.05)]
        r=convergence(levels)
        self.assertAlmostEqual(r["observed_order"],2)
        self.assertAlmostEqual(r["extrapolated"],2)
        levels[1]["value"]=1.99
        self.assertEqual(convergence(levels)["status"],"NO_MONOTONE_ORDER")

    def test_all_cli_paths_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); examples=root/"references"
            def run(*args):
                result=subprocess.run([sys.executable,"-m","ice_cfd_checks",*map(str,args)],
                    cwd=ROOT,capture_output=True,text=True)
                self.assertEqual(result.returncode,0,result.stderr)
            run("reference","--output",examples)
            cases=[("score-jump","jump_exact_example.json",[]),
                ("score-rest","hydrostatic_exact_example.csv",["--surface","1","--depth-scale","1"]),
                ("score-dam","dam_break_wet_exact_example.csv",["--left-depth","1","--right-depth","0.1","--time","0.3"]),
                ("score-thermal","thermal_exact_example.csv",[]),
                ("integrate-faces","face_flux_example.csv",[]),
                ("convergence","convergence_exact_example.json",[]),
                ("linear-response","linear_response_manufactured.json",[])]
            for name,filename,extra in cases:
                run(name,"--input",examples/filename,"--output",root/(name+".json"),*extra)
            result=subprocess.run([sys.executable,"-m","ice_cfd_checks","reference","--output",str(examples)],
                                  cwd=ROOT,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)


if __name__=="__main__":
    unittest.main()
