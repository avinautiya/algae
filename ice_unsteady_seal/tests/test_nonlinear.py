import unittest
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp

from ice_unsteady_seal.model import BaseConfig,HydraulicMap,build_base
from ice_unsteady_seal.nonlinear import NonlinearSealModel


class NonlinearVerification(unittest.TestCase):
    def test_base_is_hydraulic_equilibrium(self):
        base=build_base(BaseConfig(points=33)); model=NonlinearSealModel(base)
        hydro=model.hydraulics(0.,np.zeros(len(base.x)))
        assert_allclose(hydro["discharge"],base.config.discharge,rtol=2e-8)
        assert_allclose(hydro["velocity"],base.velocity,rtol=3e-8,atol=1e-10)
        self.assertLess(hydro["boundary_residual"],1e-8)

    def test_nonlinear_head_derivative_matches_refined_operator(self):
        base=build_base(BaseConfig(points=257)); model=NonlinearSealModel(base,hydraulic_tolerance=2e-9)
        mapping=HydraulicMap(base.x,base.velocity,1.,.5,.5,.3)
        zeros=np.zeros(len(base.x)); eps=1e-5
        derivative=(model.hydraulics(eps,zeros)["discharge"]-
                    model.hydraulics(-eps,zeros)["discharge"])/(2*eps)
        assert_allclose(mapping.direct,derivative,rtol=.003)

    def test_nonlinear_propagation_and_surrogate_ledger(self):
        base=build_base(BaseConfig(points=17)); model=NonlinearSealModel(base,hydraulic_tolerance=2e-8)
        initial=np.zeros(len(base.x)); initial[0]=1e-4
        duration=.04; inflow=1e-4
        result=model.forecast([(duration,inflow)],initial,maximum_step=.02)
        independent=solve_ivp(lambda t,y:model.rhs(y,inflow)[0],(0,duration),initial,
                              rtol=3e-7,atol=2e-11,first_step=.01,max_step=.02)
        last=result["samples"][-1]
        actual=np.concatenate(([last["lake_head_perturbation"]],result["final_bed_perturbation"][1:]))
        assert_allclose(actual,independent.y[:,-1],rtol=5e-5,atol=3e-10)
        self.assertLess(result["maximum_total_head_ledger_residual"],2e-16)


if __name__=="__main__": unittest.main()
