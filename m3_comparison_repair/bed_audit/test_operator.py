import unittest
import numpy as np
from .operator import lower_ice


class Tests(unittest.TestCase):
    def test_zero_time_preserves_every_cell_even_below_reservoir_floor(self):
        z=np.array([[-2.,-.99],[-3.,0.]])
        mask=np.array([[False,True],[False,True]])
        r=lower_ice(z,mask,np.ones_like(z)*100,0.,.01)
        np.testing.assert_array_equal(z,r['bed'])
        self.assertEqual(r['ice_volume_m3'],0.)

    def test_latent_geometry_and_nonice(self):
        z=np.array([[-2.,-.99],[-3.,0.]])
        mask=np.array([[False,True],[False,True]])
        r=lower_ice(z,mask,np.ones_like(z)*100,21600.,.01)
        self.assertTrue(r['nonice_unchanged'])
        self.assertAlmostEqual(r['geometric_latent_J'],r['heat_supplied_J'],places=6)
        self.assertAlmostEqual(r['ice_volume_m3'],r['requested_volume_m3'],places=14)

    def test_floor_is_explicit_refusal_not_heat_erasure(self):
        with self.assertRaisesRegex(RuntimeError,'floor'):
            lower_ice([[-.999]],[[True]],[[100.]],21600.,1.,minimum_ice_elevation_m=-1.)

    def test_invalid(self):
        for dt in (-1.,True,float('nan')):
            with self.assertRaises(ValueError):lower_ice([[0.]],[[True]],[[1.]],dt,1.)
        with self.assertRaises(ValueError):lower_ice([[0.]],[[1]],[[1.]],1.,1.)
        with self.assertRaises(ValueError):lower_ice([[0.]],[[True]],[[-1.]],1.,1.)


if __name__=='__main__':unittest.main()
