import unittest
from cliff_geometry import signed_area, convex_polygon, retreat_by_energy, segment_exposure

class GeometryTests(unittest.TestCase):
    box=((0.,0.),(4.,0.),(4.,3.),(0.,3.))
    def test_signed_area(self):
        self.assertEqual(signed_area(self.box),12.)
        self.assertEqual(signed_area(self.box[::-1]),-12.)
    def test_map_coordinate_area(self):
        self.assertEqual(signed_area(tuple((x+500000,z+4000000) for x,z in self.box)),12.)
    def test_normal_face_retreat(self):
        r=retreat_by_energy(self.box,(1,0),3.,rho_ice=1,latent_heat=1)
        self.assertAlmostEqual(r.normal_retreat_m,1.)
        self.assertAlmostEqual(r.removed_ice_m3,3.)
        self.assertAlmostEqual(r.energy_residual_J,0.)
    def test_floor_lowering_is_distinct(self):
        r=retreat_by_energy(self.box,(0,-1),3.,rho_ice=1,latent_heat=1)
        self.assertAlmostEqual(r.normal_retreat_m,.75)
        self.assertAlmostEqual(max(z for x,z in r.vertices),2.25)
    def test_triangle_cut(self):
        # Similar-triangle remnant: x>=1 leaves triangle of area 2.
        r=retreat_by_energy(((0,0),(3,0),(0,3)),(1,0),2.5,rho_ice=1,latent_heat=1)
        self.assertAlmostEqual(r.normal_retreat_m,1.)
        self.assertAlmostEqual(signed_area(r.vertices),2.)
    def test_width_scaling(self):
        r=retreat_by_energy(self.box,(1,0),6.,width_m=2,rho_ice=1,latent_heat=1)
        self.assertAlmostEqual(r.normal_retreat_m,1.)
    def test_zero(self):
        self.assertEqual(retreat_by_energy(self.box,(1,0),0).vertices,self.box)
    def test_exposure(self):
        r=segment_exposure((0,0),(0,10),6,2)
        self.assertEqual(r,dict(debris_covered_m=2.,submerged_bare_m=4.,exposed_bare_m=4.,total_m=10.))
    def test_wet_debris_no_doublecount(self):
        r=segment_exposure((0,0),(3,4),1,2)
        self.assertEqual(r['submerged_bare_m'],0.)
        self.assertAlmostEqual(sum(r[k] for k in ('debris_covered_m','submerged_bare_m','exposed_bare_m')),5.)
    def test_horizontal_segment(self):
        self.assertEqual(segment_exposure((0,1),(3,1),2,0)['submerged_bare_m'],3.)
    def test_exhaustion_refuses(self):
        with self.assertRaises(ValueError): retreat_by_energy(self.box,(1,0),12.,rho_ice=1,latent_heat=1)
    def test_invalid(self):
        for p in (((0,0),(1,1),(0,1),(1,0)),((0,0),(2,0),(1,.2),(2,2),(0,2))):
            with self.assertRaises(ValueError): convex_polygon(p)
        for n in ((0,0),(float('nan'),1)):
            with self.assertRaises(ValueError): retreat_by_energy(self.box,n,1)

if __name__=='__main__': unittest.main()
