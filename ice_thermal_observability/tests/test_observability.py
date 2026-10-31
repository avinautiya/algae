import copy
import itertools
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.optimize import brentq
from ice_thermal_observability.channel import Channel,CrossStream,plug_cell_bulk_kernel,modal_witness
from ice_thermal_observability.recovery import sharp_interval
from ice_thermal_observability.forecast import demo_case,predict,write_forecast


def vertices(c,e,d,eps,lo,hi):
    """Independent exhaustive intersection of n active faces; tiny cases only."""
    n=len(c)
    a=np.vstack([np.eye(n),-np.eye(n),e,-e])
    b=np.r_[hi,-lo,d+eps,-d+eps]
    values=[]
    for ids in itertools.combinations(range(len(b)),n):
        aa=a[list(ids)];bb=b[list(ids)]
        if abs(np.linalg.det(aa))<1e-12:
            continue
        x=np.linalg.solve(aa,bb)
        if np.max(a@x-b)<1e-9:
            values.append(float(c@x))
    if not values:
        raise ArithmeticError("No independently enumerated vertex")
    return min(values),max(values)


class ThermalPhysicsTests(unittest.TestCase):
    def test_exact_mean_flow_and_capacity(self):
        for velocity in ("plug","shear_free"):
            for stretch in (1,1.5,2):
                m=CrossStream(Channel(velocity=velocity),33,stretch)
                self.assertAlmostEqual(float(m.u@m.dy),m.channel.speed_m_s*m.channel.depth_m,places=16)
                self.assertAlmostEqual(float(m.mass.sum()*m.channel.width_m),m.channel.capacity_W_K,places=12)

    def test_semigroup_positivity_and_composition(self):
        m=CrossStream(Channel(),40)
        f=np.linspace(0,.2,40)
        direct=m.propagate(f,.07)
        composed=m.propagate(m.propagate(f,.03),.04)
        np.testing.assert_allclose(direct,composed,rtol=1e-11,atol=1e-13)
        self.assertGreaterEqual(float(direct.min()),-1e-12)
        self.assertLessEqual(float(direct.max()),.2+1e-12)
        self.assertLess(m.bulk(f,.07),m.bulk(f,.03))

    def test_wall_flux_integral_matches_regional_heat(self):
        m=CrossStream(Channel(),48)
        f=.1*np.sin(np.pi*m.y/(2*m.channel.depth_m))
        integrated=quad(lambda x:m.wall_flux_W_m2(f,x)*m.channel.width_m,.04,.12,
                        epsabs=1e-11,epsrel=1e-10)[0]
        self.assertAlmostEqual(integrated,float(m.heat_kernel_W_K(.04,.12)@f),places=10)

    def test_independent_analytical_kernel_refinement(self):
        errors=[]
        for n in (16,32,64,128):
            m=CrossStream(Channel(velocity="plug"),n)
            f=.08+.04*np.cos(np.pi*m.y/m.channel.depth_m)
            analytic=plug_cell_bulk_kernel(m.channel,m.edges,.04,256)
            more=plug_cell_bulk_kernel(m.channel,m.edges,.04,512)
            np.testing.assert_allclose(analytic,more,rtol=0,atol=1e-15)
            errors.append(abs(m.bulk(f,.04)-float(analytic@f)))
        self.assertTrue(all(a/b>3.5 for a,b in zip(errors[:-1],errors[1:])),errors)

    def test_no_slip_spectrum_independent_shooting(self):
        ch=Channel();fine=CrossStream(ch,128)
        scale=ch.speed_m_s*ch.depth_m**2/ch.diffusivity
        def boundary(value):
            solution=solve_ivp(lambda y,z:[z[1],-value*1.5*(2*y-y*y)*z[0]],
                               [0,1],[0.,1.],rtol=2e-11,atol=1e-12)
            return float(solution.y[1,-1])
        exact=[]
        for guess in fine.eigenvalues[:2]*scale:
            exact.append(brentq(boundary,.8*guess,1.2*guess,xtol=1e-11))
        errors=[]
        for n in (32,64,128):
            m=CrossStream(ch,n)
            errors.append(float(np.max(np.abs(m.eigenvalues[:2]*scale/np.array(exact)-1))))
        self.assertTrue(errors[0]/errors[1]>3.5 and errors[1]/errors[2]>3.5,errors)

    def test_all_thermal_modes_couple_to_bulk_discrete_identity(self):
        m=CrossStream(Channel(),48)
        for i in range(12):
            phi=m.eigenvectors[:,i]/m.root_mass
            bulk_numerator=float(m.mass@phi)
            boundary=m.bed_conductance*phi[0]/m.eigenvalues[i]
            self.assertNotEqual(bulk_numerator,0)
            self.assertAlmostEqual(bulk_numerator,boundary,places=10)

    def test_source_response_independent_ode_and_ledger(self):
        m=CrossStream(Channel(),16)
        f=np.full(16,.08);q=120.0;end=.08
        def rhs(x,state):
            t=state[:-1]
            kt=m.diagonal_K*t
            kt[:-1]+=m.offdiagonal_K*t[1:]
            kt[1:]+=m.offdiagonal_K*t[:-1]
            derivative=-kt/m.mass
            derivative[-1]+=q/m.mass[-1]
            return np.r_[derivative,m.bed_conductance*t[0]]
        solution=solve_ivp(rhs,[0,end],np.r_[f,0],method="BDF",rtol=1e-10,atol=1e-12)
        self.assertTrue(solution.success)
        direct=m.propagate(f,end)+q*m.surface_response(0,end,end)
        np.testing.assert_allclose(direct,solution.y[:-1,-1],rtol=2e-8,atol=1e-10)
        supply=q*end
        heat=float(m.mass@(f-direct)+supply)
        self.assertAlmostEqual(heat,float(solution.y[-1,-1]),places=8)

    def test_source_coefficients_causal_and_nonnegative(self):
        m=CrossStream(Channel(),40)
        self.assertEqual(m.surface_heat_coefficient(.15,.2,.04,.12),0)
        self.assertGreaterEqual(m.surface_heat_coefficient(.02,.03,.04,.12),-1e-12)
        self.assertGreater(m.surface_heat_coefficient(.04,.12,.04,.12),0)
        self.assertLess(m.surface_heat_coefficient(.04,.12,.04,.12),m.channel.width_m*.08)

    def test_smooth_continuum_counterexample(self):
        ch=Channel(velocity="plug")
        w=modal_witness(ch,[.02,.04],.06,.12,.2)
        self.assertLess(max(abs(v) for v in w["calibration_difference_K"]),1e-14)
        self.assertGreater(w["regional_power_difference_W"],1e-6)
        y=np.linspace(0,ch.depth_m,10001)
        basis=np.sin(np.outer(y,(np.arange(4)+.5)*np.pi/ch.depth_m))
        for name in ("minus_coefficients_K","plus_coefficients_K"):
            t=basis@np.array(w[name])
            self.assertGreaterEqual(float(t.min()),-1e-14)
            self.assertLessEqual(float(t.max()),.2+1e-14)


class RecoveryTests(unittest.TestCase):
    def test_independent_vertex_comparison_exact_and_noisy(self):
        rng=np.random.default_rng(717)
        for n,m in ((3,1),(4,2)):
            for noisy in (False,True):
                for _ in range(8):
                    e=rng.uniform(0,1,(m,n));e/=e.sum(axis=1)[:,None]
                    true=rng.uniform(.05,.15,n)
                    d=e@true;eps=np.full(m,.002 if noisy else 0)
                    lo=np.zeros(n);hi=np.full(n,.2);c=rng.normal(size=n)
                    r=sharp_interval(c,e,d,eps,lo,hi)
                    lower,upper=vertices(c,e,d,eps,lo,hi)
                    self.assertAlmostEqual(r["lower"],lower,places=9)
                    self.assertAlmostEqual(r["upper"],upper,places=9)
                    self.assertLess(abs(r["upper_certificate"]["dual_gap"]),1e-10)

    def test_endpoint_information_collapses_range(self):
        m=CrossStream(Channel(),32);f=np.full(32,.1)
        e=np.array([m.bulk_kernel(x) for x in (.04,.12)])
        r=sharp_interval(m.heat_kernel_W_K(.04,.12),e,e@f,np.zeros(2),np.zeros(32),np.full(32,.2))
        self.assertLess(r["upper"]-r["lower"],1e-9)

    def test_upstream_information_leaves_attainable_range(self):
        m=CrossStream(Channel(),32);f=np.full(32,.1)
        e=np.array([m.bulk_kernel(x) for x in (.02,.04)])
        c=m.heat_kernel_W_K(.06,.12)
        r=sharp_interval(c,e,e@f,np.zeros(2),np.zeros(32),np.full(32,.2))
        self.assertGreater(r["upper"]-r["lower"],.01)
        for name in ("lower_witness_variables","upper_witness_variables"):
            np.testing.assert_allclose(e@np.array(r[name]),e@f,rtol=0,atol=1e-9)

    def test_more_information_and_less_noise_do_not_widen(self):
        m=CrossStream(Channel(),32);f=np.full(32,.1);c=m.heat_kernel_W_K(.06,.12)
        e=np.array([m.bulk_kernel(x) for x in (.02,.04)])
        widths=[]
        for rows,noise in ((1,.001),(2,.001),(2,0)):
            r=sharp_interval(c,e[:rows],e[:rows]@f,np.full(rows,noise),np.zeros(32),np.full(32,.2))
            widths.append(r["upper"]-r["lower"])
        self.assertTrue(widths[0]>=widths[1]>=widths[2],widths)

    def test_target_weighted_profile_moment_collapses_range(self):
        m=CrossStream(Channel(),32);f=np.full(32,.1);c=m.heat_kernel_W_K(.06,.12)
        e=np.array([c/c.sum()])
        r=sharp_interval(c,e,e@f,np.zeros(1),np.zeros(32),np.full(32,.2))
        self.assertLess(r["upper"]-r["lower"],1e-9)

    def test_infeasible_measurement_rejected(self):
        with self.assertRaisesRegex(ValueError,"infeasible"):
            sharp_interval(np.ones(4),np.full((1,4),.25),np.array([2.]),np.zeros(1),np.zeros(4),np.ones(4))

    def test_midpoint_is_minimax_on_attaining_pair(self):
        c=np.array([1.,2.,3.]);e=np.full((1,3),1/3)
        r=sharp_interval(c,e,np.array([.1]),np.zeros(1),np.zeros(3),np.full(3,.2))
        for offset in np.linspace(-.2,.2,21):
            z=r["minimax"]+offset
            self.assertGreaterEqual(max(abs(z-r["lower"]),abs(z-r["upper"]))+1e-14,r["worst_case_absolute_error"])


class ForecastTests(unittest.TestCase):
    def test_demo_energy_replay_and_scope(self):
        result=predict(demo_case())
        self.assertFalse(result["scope_flags"])
        region=result["regions"][0]
        for witness in ("lower_witness","upper_witness"):
            w=region[witness]
            self.assertAlmostEqual(w["inlet_region_sensible_W"]-w["outlet_region_sensible_W"]+w["supplied_surface_W"],w["bed_heat_W"],places=12)
            self.assertGreaterEqual(w["min_reported_temperature_K"],-1e-10)

    def test_surface_uncertainty_and_source_energy(self):
        case=demo_case()
        case["surface_heat_segments"]=[dict(start_m=.04,stop_m=.12,lower_W_m2=0.,upper_W_m2=100.,provenance="Synthetic heating interval")]
        result=predict(case)["regions"][0]
        base=predict(demo_case())["regions"][0]
        self.assertAlmostEqual(result["power_interval"]["lower"],base["power_interval"]["lower"],places=8)
        self.assertGreater(result["power_interval"]["upper"],base["power_interval"]["upper"])
        w=result["upper_witness"]
        self.assertGreater(w["supplied_surface_W"],0)
        self.assertGreaterEqual(w["bed_heat_W"],0)

    def test_future_observations_rejected(self):
        case=demo_case();case["observations"][0]["x_m"]=.05
        with self.assertRaisesRegex(ValueError,"Future"):
            predict(case)

    def test_unknown_future_data_fields_rejected(self):
        case=demo_case();case["future_outlet_temperature_K"]=.1
        with self.assertRaisesRegex(ValueError,"Unknown"):
            predict(case)

    def test_overwrite_preserved_and_cli_executes(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td);case_path=path/"case.json";out=path/"forecast"
            subprocess.run([sys.executable,"-m","ice_thermal_observability","demo","--output",str(case_path)],check=True,capture_output=True)
            subprocess.run([sys.executable,"-m","ice_thermal_observability","predict","--case",str(case_path),"--output",str(out)],check=True,capture_output=True)
            data=json.loads((out/"forecast.json").read_text())
            self.assertEqual(len(data["regions"]),1)
            before=(out/"forecast.json").read_bytes()
            with self.assertRaises(FileExistsError):
                write_forecast(demo_case(),out)
            self.assertEqual((out/"forecast.json").read_bytes(),before)

    def test_geometry_screen_and_invalid_inputs(self):
        case=demo_case();case["horizon_s"]=60000
        self.assertTrue(predict(case)["scope_flags"])
        for field,value in (("initial_upper_K",-1),("horizon_s",0)):
            case=demo_case();case[field]=value
            with self.assertRaises(ValueError):predict(case)

    def test_temperature_scale_linearity(self):
        case=demo_case();base=predict(case)["regions"][0]["power_interval"]
        case["initial_upper_K"]*=2
        for o in case["observations"]:o["bulk_excess_K"]*=2
        scaled=predict(case)["regions"][0]["power_interval"]
        for key in ("lower","upper","minimax"):
            self.assertAlmostEqual(scaled[key],2*base[key],places=9)


if __name__=="__main__":
    unittest.main()
