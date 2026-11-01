import unittest
from model import Case
from topology import run

class TopologyTests(unittest.TestCase):
    def test_equal_area_and_volume_layouts(self):
        c=Case(discharge=.07)
        a=run(c,'corridor',days=.1);b=run(c,'spanning_strip',days=.1)
        self.assertAlmostEqual(a['covered_area_m2'],b['covered_area_m2'])
        self.assertAlmostEqual(a['solid_patch_volume_m3'],b['solid_patch_volume_m3'])
    def test_bare_corridor_order_and_contact(self):
        c=Case(discharge=.07)
        r=run(c,'corridor',days=12,dt=450)
        self.assertEqual(r['wet_days'],12)
        self.assertGreater(r['minimum_bare_flow_m3_s'],0)
        for state in r['history']:
            self.assertLessEqual(state['ice_bed'][0],state['ice_bed'][2]+1e-12)
            self.assertLessEqual(state['ice_bed'][0],state['ice_bed'][1]+1e-12)
    def test_insulated_limit_has_no_ice_loss(self):
        c=Case(discharge=.07,h_water=0)
        r=run(c,'corridor',days=.1)
        self.assertEqual(r['toe_melt_m3'],0.)
    def test_energy_balance(self):
        c=Case(discharge=.07)
        for layout in ('corridor','spanning_strip'):
            r=run(c,layout,days=.1)
            self.assertLess(abs(r['energy_relative_residual']),1e-12)
if __name__=='__main__':unittest.main()
