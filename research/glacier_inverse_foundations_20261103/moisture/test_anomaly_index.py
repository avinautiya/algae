import copy,unittest
import numpy as np
from anomaly_index import positive_area,prediction_ratios,development_fit
def fixture():
    return dict(theta=np.array([.05,.06,.04,.03,.08,.09]),theta_reference=.055,air=np.array([1.,3.,2.,1.,2.,4.]),
        seconds=np.ones(6),days=np.arange(6.),buried=np.ones(6),gradient=np.ones(6),slices=[slice(0,2),slice(2,4),slice(4,6)])
class AnomalyTests(unittest.TestCase):
    def test_exact_positive_linear_area(self):
        self.assertEqual(positive_area(-2,2,600),300);self.assertEqual(positive_area(1,3,600),1200)
    def test_constant_offset_invariance_entire_calibrated_prediction(self):
        a=fixture();b=copy.deepcopy(a);b['theta']+=.031;b['theta_reference']+=.031
        np.testing.assert_allclose(prediction_ratios(a,'anomaly_exp_air',5),prediction_ratios(b,'anomaly_exp_air',5),rtol=1e-14)
    def test_reference_independence(self):
        a=fixture();b=copy.deepcopy(a);b['theta_reference']=-10.
        np.testing.assert_allclose(prediction_ratios(a,'anomaly_exp_air',5),prediction_ratios(b,'anomaly_exp_air',5),rtol=1e-14)
    def test_centered_alternative_also_offset_invariant(self):
        a=fixture();b=copy.deepcopy(a);b['theta']+=.031;b['theta_reference']+=.031
        np.testing.assert_allclose(prediction_ratios(a,'centered_softplus_air',10),prediction_ratios(b,'centered_softplus_air',10),rtol=1e-14)
    def test_ratio_model_not_offset_invariant(self):
        a=fixture();b=copy.deepcopy(a);b['theta']+=.031;b['theta_reference']+=.031
        self.assertGreater(np.max(abs(prediction_ratios(a,'shifted_ratio_air',1)-prediction_ratios(b,'shifted_ratio_air',1))),.01)
    def test_development_selection_also_offset_invariant(self):
        a=fixture();b=copy.deepcopy(a);b['theta']+=.031;b['theta_reference']+=.031
        y=10*prediction_ratios(a,'anomaly_exp_air',5)
        aa=development_fit(a,10,y,'anomaly_exp_air',[0,5,10]);bb=development_fit(b,10,y,'anomaly_exp_air',[0,5,10])
        self.assertEqual(aa['selected']['parameter'],bb['selected']['parameter'])
    def test_bounded_drift_factor(self):
        a=fixture();b=copy.deepcopy(a);b['theta']+=np.array([-.01,-.01,.01,.01,-.01,.01])
        ratio=prediction_ratios(b,'anomaly_exp_air',5)/prediction_ratios(a,'anomaly_exp_air',5)
        self.assertTrue(np.all(ratio<=np.exp(.1)+1e-14));self.assertTrue(np.all(ratio>=np.exp(-.1)-1e-14))
if __name__=='__main__':unittest.main()
