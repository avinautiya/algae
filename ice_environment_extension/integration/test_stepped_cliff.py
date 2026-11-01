import unittest
from dataclasses import replace
from stepped_cliff import SteppedCliff, advance_latent, face_partition

def shoelace(p):
    return abs(sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(p,p[1:]+p[:1])))/2

class StepTests(unittest.TestCase):
    def setUp(self):self.c=SteppedCliff(-2.,5.,-4.,0.,0.,3.,2.)

    def test_nonconvex_area_independently(self):
        self.assertAlmostEqual(shoelace(self.c.polygon())*2,self.c.ice_volume_m3)

    def test_two_energy_channels_disjoint(self):
        e=917*334000
        new,r=advance_latent(self.c,e*.6,e*.3)
        self.assertAlmostEqual(new.face_x_m,.1)
        self.assertAlmostEqual(new.floor_z_m,-.3/4.2)
        self.assertAlmostEqual(shoelace(self.c.polygon())*2-shoelace(new.polygon())*2,.9)
        self.assertAlmostEqual(r['volume_residual_m3'],0.,places=12)

    def test_floor_is_not_face_retreat(self):
        new,r=advance_latent(self.c,0.,1e6)
        self.assertEqual(new.face_x_m,self.c.face_x_m)
        self.assertLess(new.floor_z_m,self.c.floor_z_m)
        face,r=advance_latent(self.c,1e6,0.)
        self.assertEqual(face.floor_z_m,self.c.floor_z_m)
        self.assertGreater(face.face_x_m,self.c.face_x_m)

    def test_exposure_partition(self):
        p=face_partition(self.c,1.5,.5)
        self.assertEqual(p['debris_area_m2'],1.)
        self.assertEqual(p['bare_water_area_m2'],2.)
        self.assertEqual(p['bare_air_area_m2'],3.)
        self.assertEqual(sum(v for k,v in p.items() if k!='total_face_area_m2'),6.)

    def test_invalid_energy_and_domain(self):
        for e in (-1.,float('nan'),True):
            with self.assertRaises(ValueError):advance_latent(self.c,e,0.)
        with self.assertRaises(ValueError):advance_latent(self.c,1e15,0.)
        with self.assertRaises(ValueError):replace(self.c,floor_z_m=4.).validate()

    def test_zero_exhaustion_and_mask_extremes(self):
        new,r=advance_latent(self.c,0.,0.)
        self.assertEqual(new,self.c)
        with self.assertRaises(ValueError):advance_latent(self.c,0.,1e15)
        self.assertEqual(face_partition(self.c,-10.,-10.)['bare_air_area_m2'],6.)
        self.assertEqual(face_partition(self.c,10.,-10.)['bare_water_area_m2'],6.)
        self.assertEqual(face_partition(self.c,10.,10.)['debris_area_m2'],6.)
        for v in (0.,float('nan')):
            with self.assertRaises(ValueError):advance_latent(self.c,0.,0.,rho_ice=v)

    def test_translation_invariance_at_map_coordinates(self):
        c=replace(self.c,upstream_edge_x_m=358941.,ice_back_edge_x_m=358948.,
                  face_x_m=358943.,ice_bottom_z_m=4076.,floor_z_m=4080.,crest_z_m=4083.)
        self.assertEqual(c.ice_volume_m3,self.c.ice_volume_m3)
        new,r=advance_latent(c,1e6,2e6)
        local=[(x-new.upstream_edge_x_m,z-new.ice_bottom_z_m) for x,z in new.polygon()]
        self.assertAlmostEqual(shoelace(local)*new.extrusion_width_m,new.ice_volume_m3,places=10)

if __name__=='__main__':unittest.main()
