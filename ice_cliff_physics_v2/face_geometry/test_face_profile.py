from dataclasses import replace
import math
import unittest
from face_profile import Band, Profile, advance, face_partition, front_at, RHO_ICE, LATENT


def polygon_area(points):
    # Independent shoelace quadrature, shifted to avoid absolute-datum cancellation.
    x0,z0=points[0]
    return .5*abs(math.fsum((x-x0)*(zn-z0)-(xn-x0)*(z-z0)
                           for (x,z),(xn,zn) in zip(points,points[1:]+points[:1])))


def base(n=4):
    return Profile.uniform(0.,10.,-2.,2.,0.,2.,1.,n)


class GeometryTests(unittest.TestCase):
    def test_nonconvex_volume_matches_independent_polygon(self):
        for fronts in [[2.,2.,2.],[4.,2.,3.],[2.,5.,1.]]:
            g=Profile(0.,10.,-2.,3.,2.,tuple(Band(i,i+1,x) for i,x in enumerate(fronts)))
            self.assertAlmostEqual(g.ice_volume_m3,polygon_area(g.polygon())*2.)

    def test_uniform_stepped_limit_for_any_band_count(self):
        for n in [1,3,17]:
            g=base(n);face_energy=2e6;floor_energy=3e6
            es=[face_energy*(b.z_high_m-b.z_low_m)/2 for b in g.bands]
            after,ledger=advance(g,es,floor_energy)
            x=2+face_energy/(RHO_ICE*LATENT*2)
            z=-floor_energy/(RHO_ICE*LATENT*x)
            self.assertTrue(all(abs(b.front_x_m-x)<1e-14 for b in after.bands))
            self.assertAlmostEqual(after.floor_z_m,z)
            expected=(10-0)*(z+2)+(10-x)*(2-z)
            self.assertAlmostEqual(after.ice_volume_m3,expected)
            self.assertLess(abs(ledger['volume_residual_m3']),1e-14)

    def test_face_and_floor_losses_are_disjoint(self):
        g=base(3);energy=[1e6,5e5,2e5];floor=3e6
        after,ledger=advance(g,energy,floor)
        expected=(sum(energy)+floor)/(RHO_ICE*LATENT)
        self.assertAlmostEqual(polygon_area(g.polygon())-polygon_area(after.polygon()),expected)
        self.assertAlmostEqual(ledger['stable_geometry_volume_loss_m3'],expected)
        self.assertAlmostEqual(ledger['floor_area_used_m2'],after.bands[0].front_x_m)

    def test_floor_growth_extends_lowest_band_without_extra_face_energy(self):
        g=base(2);after,ledger=advance(g,[0.,0.],RHO_ICE*LATENT*.2)
        self.assertEqual(after.bands[0].front_x_m,g.bands[0].front_x_m)
        self.assertAlmostEqual(after.vertical_face_area_m2-g.vertical_face_area_m2,.1)
        self.assertAlmostEqual(after.ice_volume_m3,g.ice_volume_m3-.2)

    def test_undercut_ceiling_not_conflated_with_vertical_face(self):
        g=Profile(0.,10.,-2.,2.,2.,(Band(0,1,4.),Band(1,2,2.)))
        self.assertEqual(g.undercut_diagnostic()['maximum_horizontal_undercut_m'],2.)
        self.assertEqual(g.undercut_diagnostic()['downward_ceiling_area_m2'],4.)
        self.assertEqual(g.vertical_face_area_m2,4.)
        self.assertEqual(g.ledges()[0]['orientation'],'downward_ceiling')

    def test_top_retreat_creates_upward_terrace(self):
        g=Profile(0.,10.,-2.,2.,1.,(Band(0,1,2.),Band(1,2,4.)))
        self.assertEqual(g.undercut_diagnostic()['maximum_horizontal_undercut_m'],0.)
        self.assertEqual(g.ledges()[0]['orientation'],'upward_terrace')

    def test_datum_translation_does_not_change_signed_melt(self):
        g=base(2)
        shifted=Profile(100.,110.,-52.,-48.,1.,tuple(Band(b.z_low_m-50,b.z_high_m-50,b.front_x_m+100) for b in g.bands))
        a,la=advance(g,[1e6,2e6],1e6);b,lb=advance(shifted,[1e6,2e6],1e6)
        self.assertAlmostEqual(a.ice_volume_m3,b.ice_volume_m3)
        self.assertAlmostEqual(la['total_melt_m3'],lb['total_melt_m3'])
        self.assertAlmostEqual(a.floor_z_m,b.floor_z_m+50.)

    def test_invalid_geometry_and_exhaustion_refused_not_clipped(self):
        with self.assertRaises(ValueError):Profile(0,10,-2,2,1,(Band(0,1,2),Band(1.01,2,2)))
        with self.assertRaises(ValueError):Profile(0,10,-2,2,1,())
        with self.assertRaises(ValueError):advance(base(1),[RHO_ICE*LATENT*20],0.)
        with self.assertRaises(ValueError):advance(base(1),[0.],RHO_ICE*LATENT*10)
        with self.assertRaises(ValueError):advance(base(1),[-1.],0.)


class PartitionTests(unittest.TestCase):
    def test_exact_subband_water_and_debris_intersections(self):
        g=base(2);p=face_partition(g,.75,[(.5,1.5)])
        self.assertEqual(p['totals'],dict(bare_water_m2=.5,bare_air_m2=.5,debris_water_m2=.25,debris_air_m2=.75))
        self.assertEqual(p['area_residual_m2'],0.)

    def test_full_and_empty_contact_and_clipped_masks(self):
        g=base(3)
        self.assertAlmostEqual(face_partition(g,10,[(-10,10)])['totals']['debris_water_m2'],2.)
        self.assertAlmostEqual(face_partition(g,-10,[])['totals']['bare_air_m2'],2.)

    def test_duplicate_cover_and_bad_stage_refused(self):
        with self.assertRaises(ValueError):face_partition(base(),1.,[(0,1),(.5,2)])
        with self.assertRaises(ValueError):face_partition(base(),float('nan'))

    def test_query_front_at_boundaries_and_outside(self):
        g=Profile(0,10,-2,2,1,(Band(0,1,3),Band(1,2,2)))
        self.assertEqual(front_at(g,0),3)
        self.assertEqual(front_at(g,1),2)
        self.assertEqual(front_at(g,2),2)
        with self.assertRaises(ValueError):front_at(g,-1)


class RefinementTests(unittest.TestCase):
    def test_nonuniform_linear_flux_profile_refines_at_first_order(self):
        # Prescribed q(z)=100+50*z, dt1000s. Exact dx(z)=q(z)*dt/(rho*L).
        errors=[]
        for n in [4,8,16]:
            g=base(n);es=[(100+25*(b.z_low_m+b.z_high_m))*1000*(b.z_high_m-b.z_low_m) for b in g.bands]
            after,ledger=advance(g,es,0.)
            self.assertAlmostEqual(ledger['face_melt_m3'],300000/(RHO_ICE*LATENT))
            error=max(abs(front_at(after,b.z_low_m+.25*(b.z_high_m-b.z_low_m))-
                          (2+(100+50*(b.z_low_m+.25*(b.z_high_m-b.z_low_m)))*1000/(RHO_ICE*LATENT))) for b in g.bands)
            errors.append(error)
        self.assertAlmostEqual(errors[0]/errors[1],2.,places=8)
        self.assertAlmostEqual(errors[1]/errors[2],2.,places=8)

    def test_moving_floor_case_approaches_independent_continuum_geometry(self):
        from run_example import analytic,run_case
        exact=analytic();errors=[]
        for dt in [3600.,1800.,900.]:
            out=run_case(16,dt)
            self.assertAlmostEqual(out['final_floor_area_m2'],exact['wet_front_x_m'],places=12)
            errors.append(abs(out['total_melt_m3']-exact['total_melt_m3']))
        self.assertTrue(1.95<errors[0]/errors[1]<2.05)
        self.assertTrue(1.95<errors[1]/errors[2]<2.05)

    def test_uniform_model_has_own_distinct_consistent_continuum_limit(self):
        from run_example import analytic,analytic_uniform,run_case
        expected=analytic_uniform();errors=[]
        self.assertNotAlmostEqual(expected['total_melt_m3'],analytic()['total_melt_m3'],places=3)
        for dt in [3600.,1800.,900.]:
            out=run_case(1,dt)
            errors.append(abs(out['total_melt_m3']-expected['total_melt_m3']))
        self.assertTrue(1.9<errors[0]/errors[1]<2.1)
        self.assertTrue(1.9<errors[1]/errors[2]<2.1)


if __name__=='__main__':unittest.main()
