import unittest
import numpy as np
from scipy.integrate import solve_ivp
from numpy.testing import assert_allclose

from ice_unsteady_seal.model import BaseConfig,build_base
from ice_unsteady_seal.matched import MatchedSealModel


class MatchedVerification(unittest.TestCase):
    def test_positive_contraction_shape_generator_and_kernel(self):
        model=MatchedSealModel(build_base(BaseConfig(points=65)))
        a=model.transport_shape
        off=a.copy(); np.fill_diagonal(off,0)
        self.assertGreaterEqual(float(np.min(off)),-1e-13)
        self.assertLessEqual(float(np.max(a.sum(axis=1))),1e-12)
        self.assertLessEqual(float(np.max(model.flux_to_shape)),1e-12)
        kernel=model.renewal_kernel([0.,.1,.5,1.])
        self.assertTrue(np.all(kernel>0))

    def test_translation_identity_converges_under_spatial_refinement(self):
        defects=[]
        for n in [65,129,257]:
            model=MatchedSealModel(build_base(BaseConfig(points=n,left=-12,right=6)))
            defects.append(abs(model.matching_diagnostics['hydraulic_translation_identity_defect']))
        self.assertLess(defects[1],defects[0]/3)
        self.assertLess(defects[2],defects[1]/3)

    def test_matched_ledger_and_independent_time_propagation(self):
        model=MatchedSealModel(build_base(BaseConfig(points=33)),2.)
        state=model.initial_state(relative_head=1e-4)
        result=model.forecast([(.3,.001)],state,sample_step=.1)
        independent=solve_ivp(lambda t,y:model.operator@y+model.input*.001,
                              (0,.3),state,rtol=2e-10,atol=1e-13)
        last=result['samples'][-1]
        actual=np.concatenate(([last['outer_seal_relative_head']],result['final_shape'][1:-1]))
        assert_allclose(actual,independent.y[:,-1],rtol=2e-8,atol=2e-12)
        self.assertLess(result['maximum_matched_limit_storage_ledger_residual'],1e-16)
        # Inner maximum and outer line-intersection height are distinct.
        self.assertIn('inner_crest_relative_head',last)
        self.assertIn('finite_crop_hydraulic_replay_difference',last)

    def test_noncanonical_case_is_rejected(self):
        base=build_base(BaseConfig(points=33,alpha=.4))
        with self.assertRaises(ValueError): MatchedSealModel(base)


if __name__=='__main__': unittest.main()
