import unittest
from dataclasses import replace
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad, solve_ivp

from ice_unsteady_seal.model import (
    BaseConfig, HydraulicMap, SealModel, build_base, hydraulic_coefficients,
    total_water_head,
)


class HydraulicVerification(unittest.TestCase):
    def test_uniform_matched_zero_shape_memory(self):
        x=np.linspace(0,2,65); v0=.7; q0=.9; alpha=.5; beta=.5; fr=.2
        mapping=HydraulicMap(x,np.full(len(x),v0),q0,alpha,beta,fr)
        ell=(1-alpha)*v0/((3-alpha)*q0)
        depth=(q0/v0)**beta
        expected=(2*beta*depth+(1-alpha)*fr**2*v0**2)/((3-alpha)*q0)
        j,e,v=mapping.solve(.002,np.zeros(len(x)))
        assert_allclose(j,.002/expected,rtol=2e-13)
        assert_allclose(v,ell*j,rtol=5e-13,atol=1e-17)
        # The upstream-reference incision must be subtracted in this gauge.
        shape_rate=-3*v0**2*(v-ell*j)
        assert_allclose(shape_rate,0,atol=1e-16)

    def test_shape_sensitivity_against_independent_quadrature(self):
        errors=[]; length=2.; v0=.7; q0=.9; alpha=.5; beta=.5; fr=.2
        d,hq,r,rv,rq=hydraulic_coefficients(v0,q0,alpha,beta,fr)
        a=float(rv/d); ell=(1-alpha)*v0/((3-alpha)*q0); z=float(hq-d*ell)
        bump=lambda x: x*x*(length-x)**2/length**4
        exact=-a/z*quad(lambda x:np.exp(-a*x)*bump(x),0,length,epsabs=1e-13)[0]
        for n in [33,65,129]:
            x=np.linspace(0,length,n)
            mapping=HydraulicMap(x,np.full(n,v0),q0,alpha,beta,fr)
            j,_,_=mapping.solve(0,bump(x))
            self.assertLess(j,0); errors.append(abs(j-exact))
        self.assertLess(errors[1],errors[0]/3.5)
        self.assertLess(errors[2],errors[1]/3.5)

    def test_critical_state_is_rejected(self):
        with self.assertRaises(ValueError):
            HydraulicMap(np.linspace(0,1,9),np.ones(9),1,.5,.5,1.)


class NonuniformVerification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base=build_base(BaseConfig(points=65))
        cls.model=SealModel(cls.base)

    def test_actual_nonuniform_seal_and_matching_diagnostics(self):
        base=self.base; c=base.config
        self.assertLess(base.crest_curvature,0)
        self.assertLess(base.velocity[0],base.crest_velocity)
        self.assertGreater(base.velocity[-1],base.crest_velocity)
        crest=np.argmin(abs(base.x))
        assert_allclose(base.bed[crest],0,atol=1e-13)
        assert_allclose(base.upstream_slope-base.crest_velocity**3/c.relative_ice_speed,0,atol=1e-13)
        _,_,r,_,_=hydraulic_coefficients(base.velocity[-1],c.discharge,c.alpha,c.beta,c.froude_coefficient)
        normal_mismatch=base.upstream_slope-base.velocity[-1]**3/c.relative_ice_speed+float(r)
        # The downstream asymptote is approximate at finite cutoff, not exact.
        self.assertLess(abs(normal_mismatch),1e-8)

    def test_same_head_crest_and_end_slopes_different_flux(self):
        x=self.base.x
        shape=np.zeros_like(x)
        inside=(x>1)&(x<3)
        z=(x[inside]-2)
        shape[inside]=1e-4*np.exp(-1/(1-z*z))
        j,_,_=self.model.hydraulics.solve(0,shape)
        self.assertLess(j,0)
        self.assertEqual(shape[np.argmin(abs(x))],0)
        self.assertTrue(np.all(shape[:3]==0) and np.all(shape[-3:]==0))

    def test_forecast_against_independent_time_integrator_and_ledger(self):
        base=build_base(BaseConfig(points=25)); model=SealModel(base,2.)
        state=model.initial_state(head=1e-4)
        forcing=1e-4; duration=.4
        result=model.forecast([(duration,forcing)],state,sample_step=.1)
        numerical=solve_ivp(lambda t,y:model.operator@y+model.input*forcing,
                            (0,duration),state,rtol=1e-10,atol=1e-13)
        exact=model.initial_state(result["samples"][-1]["lake_head_perturbation"],
                                  result["final_bed_perturbation"])
        assert_allclose(exact,numerical.y[:,-1],rtol=3e-8,atol=2e-12)
        self.assertLess(result["maximum_total_head_ledger_residual"],2e-16)
        self.assertIn("maximum_surface_stage_storage_residual",result)

    def test_nonuniform_linear_map_against_nonlinear_shooting(self):
        # Independently integrate the nonlinear spatial velocity ODE while
        # perturbing the bed. Compare its implicit head/flux derivative.
        base=build_base(BaseConfig(points=513)); c=base.config
        def bump(x):
            if abs(x-2)>=1: return 0.
            return np.exp(-1/(1-(x-2)**2))
        def bump_derivative(x):
            z=x-2
            if abs(z)>=1: return 0.
            return -2*z*np.exp(-1/(1-z*z))/(1-z*z)**2
        def shoot(q,epsilon):
            exponent=(1-c.alpha)/(3-c.alpha)
            end=base.velocity[-1]*(q/c.discharge)**exponent
            def rhs(x,y):
                vb=base.evaluate([x])[0,0]
                slope=base.upstream_slope-vb**3/c.relative_ice_speed+epsilon*bump_derivative(x)
                d,_,r,_,_=hydraulic_coefficients(y[0],q,c.alpha,c.beta,c.froude_coefficient)
                return [(slope+float(r))/float(d)]
            run=solve_ivp(rhs,(c.right,c.left),[end],rtol=2e-11,atol=2e-13,
                          max_step=.06,first_step=.01)
            self.assertTrue(run.success)
            return float(base.bed[0]+total_water_head(run.y[0,-1],q,c.beta,c.froude_coefficient))
        eps=1e-4
        dh_dq=(shoot(c.discharge+eps,0)-shoot(c.discharge-eps,0))/(2*eps)
        dh_de=(shoot(c.discharge,eps)-shoot(c.discharge,-eps))/(2*eps)
        independent=-dh_de/dh_dq
        errors=[]
        for points in [129,257,513]:
            grid=np.linspace(c.left,c.right,points)
            velocity=base.evaluate(grid)[0]
            mapping=HydraulicMap(grid,velocity,c.discharge,c.alpha,c.beta,c.froude_coefficient)
            linear=mapping.solve(0,np.asarray([bump(x) for x in grid]))[0]
            errors.append(abs(linear-independent))
        self.assertLess(errors[1],errors[0]/3)
        self.assertLess(errors[2],errors[1]/3)
        assert_allclose(linear,independent,rtol=.006,atol=1e-8)

    def test_initial_state_and_time_rules(self):
        with self.assertRaises(ValueError):
            self.model.initial_state(bed=np.ones(len(self.base.x)))
        with self.assertRaises(ValueError): self.model.forecast([(-1,0)])
        with self.assertRaises(ValueError): self.model.forecast([(1,-2)])

    def test_small_amplitude_does_not_override_critical_margin(self):
        base=build_base(BaseConfig(points=129,froude_coefficient=np.sqrt(.48)))
        model=SealModel(base,subcritical_margin=.02)
        bed=np.zeros(len(base.x))
        bed[np.argmin(abs(base.x-5))]=.015
        result=model.forecast([],model.initial_state(bed=bed))
        sample=result["samples"][0]
        self.assertLess(max(sample["relative_perturbations"].values()),.1)
        self.assertFalse(sample["inside_declared_subcritical_margin"])
        self.assertFalse(sample["within_declared_linear_screen"])


if __name__=="__main__": unittest.main()
