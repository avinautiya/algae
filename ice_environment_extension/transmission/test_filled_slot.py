import unittest
from filled_slot import filled_slot_discharge

class SlotTests(unittest.TestCase):
    def setUp(self):
        self.c=dict(upstream_stage_m=.2,downstream_stage_m=.1,floor_z_m=0.,
                    aperture_m=1e-5,width_m=.5,length_m=.1,viscosity_pa_s=.00179,
                    rho_water_kg_m3=1000.,entrance_discharge_coefficient=.6)

    def test_poiseuille_limit_and_head_ledger(self):
        r=filled_slot_discharge(**self.c)
        expected=.5*(1e-5)**3*1000*9.80665*.1/(12*.00179*.1)
        self.assertAlmostEqual(r['flow_m3_s']/expected,1.,places=8)
        self.assertAlmostEqual(r['head_residual_m'],0.,places=14)
        self.assertLess(r['reynolds_Dh'],1.)

    def test_cubic_small_gap_and_path_length(self):
        r=filled_slot_discharge(**self.c)['flow_m3_s']
        wider=filled_slot_discharge(**{**self.c,'aperture_m':2e-5})['flow_m3_s']
        longer=filled_slot_discharge(**{**self.c,'length_m':.2})['flow_m3_s']
        self.assertAlmostEqual(wider/r,8.,places=6)
        self.assertAlmostEqual(longer/r,.5,places=8)

    def test_invalid_partial_fill_and_zero_head(self):
        self.assertEqual(filled_slot_discharge(**{**self.c,'upstream_stage_m':.1})['flow_m3_s'],0.)
        with self.assertRaises(ValueError):filled_slot_discharge(**{**self.c,'downstream_stage_m':0.})
        with self.assertRaises(ValueError):filled_slot_discharge(**{**self.c,'aperture_m':.2})

if __name__=='__main__':unittest.main()
