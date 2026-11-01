import unittest
from ice_cliff_survival import coupling, UnsupportedRun
class TestUnfinished(unittest.TestCase):
    def test_mobile_grain_coupling_fails_explicitly(self):
        with self.assertRaises(UnsupportedRun): coupling.MobileGrainCoupling()
    def test_static_shield_rows(self):
        h = coupling.static_shield_rows(20, 10.0, 3.0, 62.0, 33.0, 0.5)
        self.assertEqual(len(h), 20); self.assertTrue(all(a <= b for a, b in zip(h, h[1:]))); self.assertEqual(h[-1], 1e30)      # bare above the pile; conductance rises toward the crest
