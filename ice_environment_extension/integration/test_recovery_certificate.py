import unittest
from dataclasses import replace
from recovery_certificate import BoundInputs, no_contact_bound

class BoundTests(unittest.TestCase):
    def setUp(self):
        self.c=BoundInputs(0.,0.,1.,.008,.3,0.,2.,0.)

    def test_stationary_bypass_barrier(self):
        r=no_contact_bound(self.c)
        self.assertAlmostEqual(r['stage_upper_m'],.04)
        self.assertAlmostEqual(r['no_contact_margin_m'],.26)
        self.assertEqual(r['verdict'],'contact_excluded_under_stated_bounds')

    def test_energy_threshold_and_cold_storage(self):
        energy=917*334000*2*.27
        hot=no_contact_bound(replace(self.c,toe_floor_positive_heat_upper_J=energy))
        self.assertEqual(hot['verdict'],'undetermined')
        cold=no_contact_bound(replace(self.c,toe_floor_positive_heat_upper_J=energy,
                                      toe_floor_cold_content_lower_J=energy))
        self.assertEqual(cold['verdict'],'contact_excluded_under_stated_bounds')

    def test_high_initial_storage_controls_bound(self):
        r=no_contact_bound(replace(self.c,initial_stage_m=.2))
        self.assertEqual(r['stage_upper_m'],.2)

    def test_uncertainty_only_weakens_exclusion(self):
        base=no_contact_bound(self.c)['no_contact_margin_m']
        for c in (replace(self.c,total_positive_reservoir_inflow_upper_m3_s=.04),
                  replace(self.c,bypass_rating_coefficient_lower_m12_s=.5),
                  replace(self.c,toe_floor_positive_heat_upper_J=1e7)):
            self.assertLess(no_contact_bound(c)['no_contact_margin_m'],base)

    def test_invalid_inputs(self):
        for c in (replace(self.c,minimum_contact_offset_above_ice_m=-.1),
                  replace(self.c,toe_floor_area_m2=0),
                  replace(self.c,initial_stage_m=float('nan')),
                  replace(self.c,total_positive_reservoir_inflow_upper_m3_s=True)):
            with self.assertRaises(ValueError):no_contact_bound(c)

if __name__=='__main__':unittest.main()
