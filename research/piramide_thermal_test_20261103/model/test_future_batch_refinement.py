import unittest
from future_batch_refinement import compare_block


class LongTrajectoryMetricsTests(unittest.TestCase):
    def test_fixed_gain_multiplier_without_reading_gain(self):
        baseline=dict(start='a',end='b',normalized_basal_heat_K_m=2.,raw_melt_mmwe=10.)
        current=dict(start='a',end='b',normalized_basal_heat_K_m=2.1,raw_melt_mmwe=10.5)
        r=compare_block(current,baseline)
        self.assertAlmostEqual(r['relative_fixed_version3_gain_prediction_change'],.05)
        self.assertAlmostEqual(r['fixed_version3_gain_prediction_multiplier'],1.05)

    def test_clock_or_zero_dose_mismatch_refused(self):
        a=dict(start='a',end='b',normalized_basal_heat_K_m=2.,raw_melt_mmwe=1.)
        with self.assertRaises(ValueError):compare_block(dict(a,end='c'),a)
        with self.assertRaises(ValueError):compare_block(a,dict(a,normalized_basal_heat_K_m=0.))


if __name__=='__main__':unittest.main()
