import math
import unittest
from heat_delivery import water_heat_rates

class HeatTests(unittest.TestCase):
    def test_combined_exchanger_and_ledger(self):
        r=water_heat_rates(.01,.5,100.,200.)
        expected=.5*math.exp(-300/(1000*4217*.01))
        self.assertAlmostEqual(r['outgoing_superheat_K'],expected)
        self.assertAlmostEqual(r['energy_residual_W'],0.)
        self.assertLess(r['floor_water_heat_W']+r['face_water_heat_W'],r['incoming_water_heat_W'])

    def test_infinite_floor_exhausts_water_heat(self):
        r=water_heat_rates(1e-6,1.,1e10,1e10)
        self.assertEqual(r['face_water_heat_W'],0.)
        self.assertEqual(r['outgoing_water_heat_W'],0.)
        self.assertGreater(r['floor_water_heat_W'],0.)

    def test_no_heat_without_water_or_exchange(self):
        self.assertEqual(water_heat_rates(0.,.5,1.,1.)['incoming_water_heat_W'],0.)
        r=water_heat_rates(.01,.5,0.,0.)
        self.assertEqual(r['incoming_water_heat_W'],r['outgoing_water_heat_W'])
        with self.assertRaises(ValueError):water_heat_rates(.1,-1.,1.,1.)

if __name__=='__main__':unittest.main()
