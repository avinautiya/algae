import math
import random
import unittest
from face_profile import Band, Profile, face_partition, front_at
from integration_adapter import initialize_profile, advance_with_partition
from mesh import split_at_elevations

GEOMETRY=dict(upstream_edge_x_m=0.,ice_back_edge_x_m=10.,ice_bottom_z_m=-2.,
              face_x_m=2.,floor_z_m=0.,crest_z_m=2.,extrusion_width_m=1.)


def polygon_area(points):
    return .5*abs(math.fsum(x*zn-xn*z for (x,z),(xn,zn) in zip(points,points[1:]+points[:1])))


class MeshTests(unittest.TestCase):
    def test_nonplanar_shape_volume_and_ledge_surfaces_unchanged(self):
        g=Profile(0.,10.,-2.,3.,2.,(Band(0,1,4),Band(1,2,2),Band(2,3,3)))
        new=split_at_elevations(g,[.2,2.6,.6,2.,.2,1.1])
        self.assertEqual(len(new.bands),7)
        self.assertAlmostEqual(new.ice_volume_m3,g.ice_volume_m3)
        self.assertAlmostEqual(polygon_area(new.polygon()),polygon_area(g.polygon()))
        self.assertEqual(new.ledges(),g.ledges())
        self.assertEqual(new.undercut_diagnostic(),g.undercut_diagnostic())
        for z in [0,.2,.8,1,1.1,1.9,2,2.6,3]:self.assertEqual(front_at(new,z),front_at(g,z))

    def test_initial_uniform_front_and_floor_preserved_exactly(self):
        g=initialize_profile(GEOMETRY,1);new=split_at_elevations(g,[.31,.73,1.21])
        self.assertEqual(new.floor_z_m,g.floor_z_m)
        self.assertEqual(new.floor_area_m2,g.floor_area_m2)
        self.assertEqual(new.ice_volume_m3,g.ice_volume_m3)
        self.assertTrue(all(b.front_x_m==2. for b in new.bands))

    def test_already_present_or_outside_splits_are_noops(self):
        g=initialize_profile(GEOMETRY,2)
        self.assertIs(split_at_elevations(g,[-100,0,1,1,2,100]),g)
        self.assertIs(split_at_elevations(g,[]),g)

    def test_cap_refuses_without_mutation_or_dropped_points(self):
        g=initialize_profile(GEOMETRY,1);before=repr(g);points=[.1,.2,.3]
        with self.assertRaisesRegex(ValueError,'exceeding cap'):
            split_at_elevations(g,points,max_bands=3)
        self.assertEqual(repr(g),before);self.assertEqual(points,[.1,.2,.3])
        self.assertEqual(len(split_at_elevations(g,points,max_bands=4).bands),4)
        with self.assertRaises(ValueError):split_at_elevations(initialize_profile(GEOMETRY,4),[],max_bands=3)

    def test_no_epsilon_merging_of_distinct_representable_interfaces(self):
        g=initialize_profile(GEOMETRY,1);a=.3;b=math.nextafter(a,math.inf)
        new=split_at_elevations(g,[a,b]);self.assertEqual(len(new.bands),3)
        self.assertEqual(new.bands[1].z_low_m,a);self.assertEqual(new.bands[1].z_high_m,b)

    def test_invalid_elevations_and_caps_refused(self):
        g=initialize_profile(GEOMETRY,1)
        for zs in [[math.nan],[math.inf],[True],'0.5',None]:
            with self.assertRaises(ValueError):split_at_elevations(g,zs)
        for cap in [0,-1,True,2.0]:
            with self.assertRaises(ValueError):split_at_elevations(g,[],max_bands=cap)

    def test_deterministic_many_splits_preserve_reference_polygon(self):
        rng=random.Random(516)
        g=Profile(0.,10.,-2.,5.,.7,tuple(Band(i,i+1,rng.uniform(1,9)) for i in range(5)))
        zs=[rng.uniform(0,5) for _ in range(30)];new=split_at_elevations(g,zs)
        self.assertLess(abs(new.ice_volume_m3-g.ice_volume_m3),1e-13)
        self.assertLess(abs(polygon_area(new.polygon())-polygon_area(g.polygon())),1e-13)

    def test_aligned_uniform_class_flux_independent_of_starting_grid(self):
        # Both grids align the sharp material and water interfaces before energy
        # allocation. They therefore represent the SAME classwise normal speed.
        stage=.73;cover=[(.31,1.21)];out=[]
        for n in [1,8,31]:
            g=split_at_elevations(initialize_profile(GEOMETRY,n),[stage,.31,1.21])
            p=face_partition(g,stage,cover);qdt={'bare_air':1e5,'bare_water':1e6,'covered_dry':2e5,'covered_water':3e5}
            area={'bare_air':'bare_air_m2','bare_water':'bare_water_m2','covered_dry':'debris_air_m2','covered_water':'debris_water_m2'}
            energies={k:v*p['totals'][area[k]] for k,v in qdt.items()}
            after,ledger,_=advance_with_partition(g,energies,4e5,stage,cover)
            out.append((after,ledger))
        ref,refledger=out[0]
        for g,ledger in out[1:]:
            self.assertAlmostEqual(g.ice_volume_m3,ref.ice_volume_m3)
            self.assertAlmostEqual(g.floor_z_m,ref.floor_z_m)
            self.assertAlmostEqual(ledger['total_melt_m3'],refledger['total_melt_m3'])
            for z in [0,.1,.5,.8,1.5,2]:self.assertAlmostEqual(front_at(g,z),front_at(ref,z))


if __name__=='__main__':unittest.main()
