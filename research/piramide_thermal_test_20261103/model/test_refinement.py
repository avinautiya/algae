import math,unittest
from refinement import metrics


class RefinementMetricsTests(unittest.TestCase):
    def test_index_units_and_both_calibration_interpretations(self):
        baseline=dict(training_normalized_basal_heat_K_m=2.,evaluation_normalized_basal_heat_K_m=6.,evaluation_training_heat_ratio=3.)
        r=metrics([1.,3.,2.,10.],2,baseline)
        self.assertEqual(r['training_normalized_basal_heat_K_m'],4.)
        self.assertEqual(r['evaluation_normalized_basal_heat_K_m'],12.)
        self.assertEqual(r['fixed_original_melt_gain_prediction_multiplier'],2.)
        self.assertEqual(r['same_train_target_recalibrated_gain_prediction_multiplier'],1.)
        self.assertEqual(r['relative_change_heat_ratio'],0.)

    def test_nonpositive_heat_refused_without_clipping(self):
        baseline=dict(training_normalized_basal_heat_K_m=1.,evaluation_normalized_basal_heat_K_m=1.,evaluation_training_heat_ratio=1.)
        with self.assertRaises(ValueError):metrics([-1.,1.],1,baseline)
        with self.assertRaises(ValueError):metrics([1.,-1.],1,baseline)


if __name__=='__main__':unittest.main()
