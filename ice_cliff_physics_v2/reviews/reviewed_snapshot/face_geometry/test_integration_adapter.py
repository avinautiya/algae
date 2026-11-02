import math
import unittest
from face_profile import RHO_ICE, LATENT, front_at
from integration_adapter import initialize_profile, advance_with_partition, routing_proxy


GEOMETRY=dict(upstream_edge_x_m=0.,ice_back_edge_x_m=10.,ice_bottom_z_m=-2.,
              face_x_m=2.,floor_z_m=0.,crest_z_m=2.,extrusion_width_m=1.)


class IntegrationAdapterTests(unittest.TestCase):
    def test_all_four_classes_and_floor_have_one_conserved_energy_budget(self):
        g=initialize_profile(GEOMETRY,8)
        energies=dict(bare_air=1e5,bare_water=2e5,covered_dry=3e5,covered_water=4e5)
        after,ledger,proxy=advance_with_partition(g,energies,5e5,.75,[(.5,1.5)])
        expected=(sum(energies.values())+5e5)/(RHO_ICE*LATENT)
        self.assertAlmostEqual(g.ice_volume_m3-after.ice_volume_m3,expected)
        self.assertAlmostEqual(ledger['total_melt_m3'],expected)
        self.assertEqual(ledger['allocation_residual_J'],0.)
        self.assertEqual(ledger['class_areas_before_m2'],dict(bare_air=.5,bare_water=.5,covered_dry=.75,covered_water=.25))
        self.assertEqual(proxy['face_x_m'],after.bands[0].front_x_m)
        self.assertEqual(proxy['floor_z_m'],after.floor_z_m)

    def test_uniform_all_class_flux_reproduces_planar_step(self):
        for n in [1,3,32]:
            g=initialize_profile(GEOMETRY,n);qdt=1e6
            energies=dict(bare_air=qdt*.5,bare_water=qdt*.5,covered_dry=qdt*.75,covered_water=qdt*.25)
            after,ledger,proxy=advance_with_partition(g,energies,3e5,.75,[(.5,1.5)])
            expected_x=2+qdt/(RHO_ICE*LATENT)
            expected_floor=-3e5/(RHO_ICE*LATENT*expected_x)
            self.assertTrue(all(abs(b.front_x_m-expected_x)<1e-14 for b in after.bands))
            self.assertAlmostEqual(after.floor_z_m,expected_floor)
            self.assertAlmostEqual(ledger['physical_minus_proxy_remaining_volume_m3'],0.)

    def test_nonuniform_class_forcing_refines_partition_shape(self):
        # Stage .73 cuts through initial bands. Band-average recession differs
        # only in that intersected band; integrated profile L1 error is O(dz).
        errors=[]
        for n in [8,32,128]:
            g=initialize_profile(GEOMETRY,n);stage=.73;wet_flux_dt=1e6;air_flux_dt=1e5
            after,ledger,_=advance_with_partition(g,dict(bare_water=wet_flux_dt*stage,
                bare_air=air_flux_dt*(2-stage),covered_dry=0.,covered_water=0.),0.,stage,[])
            error=0.;dz=2/10000
            for j in range(10000):
                z=(j+.5)*dz;exact=2+(wet_flux_dt if z<stage else air_flux_dt)/(RHO_ICE*LATENT)
                error+=abs(front_at(after,z)-exact)*dz
            errors.append(error)
            self.assertLess(abs(ledger['allocation_residual_J']),1e-8)
        self.assertGreater(errors[0],errors[1]);self.assertGreater(errors[1],errors[2])

    def test_proxy_volume_is_not_substituted_for_actual_profile(self):
        g=initialize_profile(GEOMETRY,4)
        after,ledger,proxy=advance_with_partition(g,dict(bare_water=1e7,bare_air=1e5,covered_dry=0.,covered_water=0.),0.,.5,[])
        self.assertNotEqual(ledger['after_profile_ice_volume_m3'],ledger['routing_proxy_ice_volume_m3'])
        self.assertEqual(ledger['after_profile_ice_volume_m3'],after.ice_volume_m3)
        self.assertFalse(ledger['routing_proxy_is_physical_volume'])

    def test_snapshot_partition_precedes_floor_area_growth(self):
        g=initialize_profile(GEOMETRY,4)
        after,ledger,_=advance_with_partition(g,dict(bare_air=0.,bare_water=0.,covered_dry=0.,covered_water=0.),1e7,0.,[])
        self.assertEqual(ledger['class_areas_before_m2']['bare_air'],2.)
        self.assertGreater(after.vertical_face_area_m2,2.)
        self.assertEqual(sum(ledger['band_face_energy_J']),0.)

    def test_positive_energy_for_missing_class_is_refused(self):
        g=initialize_profile(GEOMETRY,4)
        with self.assertRaisesRegex(ValueError,'covered_water'):
            advance_with_partition(g,dict(bare_air=0.,bare_water=0.,covered_dry=0.,covered_water=1.),0.,1.,[])

    def test_negative_nonfinite_or_unknown_class_energy_is_refused(self):
        g=initialize_profile(GEOMETRY,4)
        for bad in [dict(bare_air=-1.,bare_water=0.,covered_dry=0.,covered_water=0.),
                    dict(bare_air=math.inf,bare_water=0.,covered_dry=0.,covered_water=0.),
                    dict(bare_air=0.,bare_water=0.,covered_dry=0.)]:
            with self.assertRaises(ValueError):advance_with_partition(g,bad,0.,1.,[])


if __name__=='__main__':unittest.main()
