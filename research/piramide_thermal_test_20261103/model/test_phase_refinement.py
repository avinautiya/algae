import unittest
from phase_refinement import retreat,compare


class PhaseRefinementMetricsTests(unittest.TestCase):
    def test_ice_mass_units_and_empirical_qualification(self):
        r=retreat(917.,debris_depth_m=.5,ice_probe_gap_m=.02,ice_cell_width_m=.1,ice_depth_m=2.,multiplier=2.)
        self.assertEqual(r['raw_ice_equivalent_retreat_m'],1.)
        self.assertEqual(r['raw_retreat_over_initial_DT5_ice_gap'],50.)
        self.assertEqual(r['empirically_scaled_diagnostic']['ice_equivalent_retreat_m'],2.)
        self.assertTrue(r['empirically_scaled_diagnostic']['not_a_conserved_physical_solve'])

    def test_both_prediction_calibration_interpretations(self):
        r=compare(4.,12.,dict(raw_training_melt_mmwe=2.,raw_evaluation_melt_mmwe=6.))
        self.assertEqual(r['fixed_original_empirical_gain_prediction_multiplier'],2.)
        self.assertEqual(r['same_train_target_recalibrated_prediction_multiplier'],1.)

    def test_invalid_retreat_and_geometry_refused(self):
        with self.assertRaises(ValueError):retreat(-1.,debris_depth_m=.5,ice_probe_gap_m=.02,ice_cell_width_m=.1,ice_depth_m=2.)
        with self.assertRaises(ValueError):retreat(1.,debris_depth_m=.5,ice_probe_gap_m=0.,ice_cell_width_m=.1,ice_depth_m=2.)
        with self.assertRaises(ValueError):retreat(1.,debris_depth_m=.5,ice_probe_gap_m=.02,ice_cell_width_m=.1,ice_depth_m=2.,multiplier=True)


if __name__=='__main__':unittest.main()
