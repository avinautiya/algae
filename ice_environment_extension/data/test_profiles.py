import unittest
from summarize_profiles import bin_profile

class ProfileTests(unittest.TestCase):
    def test_filter_median_and_negative_bin(self):
        points=[(-.1,2.,0.),(-.2,4.,.1),(.1,10.,0.),(.2,12.,0.),(.3,100.,2.)]
        a=bin_profile(points,min_points=2)
        self.assertEqual([v['bin_index'] for v in a],[-1,0])
        self.assertEqual([v['median_elevation_m'] for v in a],[3.,11.])
        self.assertEqual([v['point_count'] for v in a],[2,2])

    def test_invalid_and_sparse(self):
        self.assertEqual(bin_profile([(0.,2.,0.)]),[])
        with self.assertRaises(ValueError):bin_profile([(0.,float('nan'),0.)])
        with self.assertRaises(ValueError):bin_profile([],bin_width_m=0)

if __name__=='__main__':unittest.main()
