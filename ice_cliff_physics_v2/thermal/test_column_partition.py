import unittest
from column_partition import ColumnPartition as P,remap


class PartitionTests(unittest.TestCase):
    def test_no_boundary_change_does_not_mix_hot_and_cold_pools(self):
        p=P(100,.3,.7,(1.,2.,3.),(6.,7.,8.))
        r=remap(p,1.,.3,(0.,0.,0.))
        self.assertEqual(r.partition.wet_profile_C,p.wet_profile_C)
        self.assertEqual(r.partition.dry_profile_C,p.dry_profile_C)
        self.assertEqual(r.energy_residual_J,0)

    def test_partial_wetting_conserves_each_cell(self):
        p=P(100,.25,.75,(0.,1.),(4.,5.))
        r=remap(p,1.,.5,(0.,0.))
        self.assertEqual(r.partition.wet_profile_C,(2.,3.))
        self.assertEqual(r.partition.dry_profile_C,(4.,5.))
        self.assertEqual(r.wetting_transfer_energy_J,112.5)
        self.assertEqual(r.cell_energy_residuals_J,(0.,0.))

    def test_wet_dry_wet_cycle_keeps_inventory(self):
        p=P(800,1.,0.,(2.,4.,6.),None)
        r=remap(p,1.,0.,(0.,0.,0.))
        self.assertIsNone(r.partition.wet_profile_C)
        self.assertEqual(r.partition.dry_profile_C,p.wet_profile_C)
        s=remap(r.partition,1.,1.,(0.,0.,0.))
        self.assertIsNone(s.partition.dry_profile_C)
        self.assertEqual(s.partition.wet_profile_C,p.wet_profile_C)
        self.assertEqual(s.partition.stored_energy_J,p.stored_energy_J)

    def test_expansion_imports_declared_hot_or_cold_profile(self):
        for incoming in [(0.,0.),(10.,20.)]:
            p=P(100,.5,.5,(1.,3.),(5.,7.))
            r=remap(p,2.,.5,incoming)
            expected=50*sum(incoming)
            self.assertEqual(r.imported_energy_J,expected)
            self.assertAlmostEqual(r.partition.stored_energy_J,p.stored_energy_J+expected)
            self.assertEqual(r.partition.wet_profile_C,p.wet_profile_C)

    def test_expansion_then_wetting_has_explicit_donor_average(self):
        p=P(100,.5,.5,(1.,1.),(3.,3.))
        r=remap(p,2.,1.,(9.,9.))
        # Dry donor: (.5*3+1*9)/1.5 = 7; half m2 transferred.
        self.assertEqual(r.partition.wet_profile_C,(4.,4.))
        self.assertEqual(r.partition.dry_profile_C,(7.,7.))
        self.assertEqual(r.imported_energy_J,900.)
        self.assertAlmostEqual(r.energy_residual_J,0)

    def test_empty_domain_can_receive_declared_material(self):
        p=P(100,0.,0.,None,None)
        r=remap(p,1.,.4,(2.,4.))
        self.assertEqual(r.partition.wet_profile_C,(2.,4.))
        self.assertEqual(r.partition.dry_profile_C,(2.,4.))
        self.assertAlmostEqual(r.partition.stored_energy_J,300)

    def test_canonical_area_prevents_one_ulp_false_contraction(self):
        from dataclasses import replace
        import math
        area=.07677247963474888
        p=P(100,0.,area,None,(2.,4.))
        first=remap(p,area,.012960718626439903,(0.,0.)).partition
        self.assertGreater(first.wet_area_m2+first.dry_area_m2,area)
        self.assertEqual(first.total_area_m2,area)
        # A thermal update must preserve the same canonical geometry reference.
        updated=replace(first,wet_profile_C=(1.,3.))
        second=remap(updated,area,.012956431833457863,(0.,0.))
        self.assertEqual(second.expanded_area_m2,0.)
        self.assertEqual(second.partition.total_area_m2,area)
        self.assertLess(abs(second.energy_residual_J),1e-12)
        with self.assertRaises(ValueError):remap(updated,math.nextafter(area,-math.inf),.01,(0.,0.))

    def test_canonical_area_cannot_hide_genuine_area_mismatch(self):
        with self.assertRaises(ValueError):P(100,.3,.7,(1.,),(2.,),area_reference_m2=.9).validate()
        with self.assertRaises(ValueError):P(100,0.,1.,None,(2.,),area_reference_m2=0.).validate()

    def test_invalid_and_contracting_areas_refused(self):
        p=P(100,.5,.5,(1.,1.),(2.,2.))
        for total,wet in [(.9,.5),(1.,1.1),(-1.,0.),(float('nan'),0.)]:
            with self.assertRaises(ValueError):remap(p,total,wet,(0.,0.))
        with self.assertRaises(ValueError):remap(p,1.,.5,(-1.,0.))
        with self.assertRaises(ValueError):remap(p,1.,.5,(0.,))
        with self.assertRaises(ValueError):P(100,0.,1.,(1.,),(0.,)).validate()


class CoupledThermalPartitionTests(unittest.TestCase):
    def test_switching_and_expansion_close_full_energy_without_global_pool_mixing(self):
        from resolved_slab import SlabState,ArealSlabParameters,advance as slab_step
        from wet_debris import ThermalState,ThermalParameters,advance as water_only
        p=P(90000.,.3,.7,tuple([1.]*8),tuple([8.]*8))
        initial=p.stored_energy_J;source=imported=ice=0.
        for total,wet in [(1.,.3),(1.,.8),(1.,0.),(1.2,1.2),(1.2,.3)]:
            rem=remap(p,total,wet,tuple([3.]*8));p=rem.partition;imported+=rem.imported_energy_J
            wet_profile=p.wet_profile_C;dry_profile=p.dry_profile_C
            if p.wet_area_m2>0:
                r=slab_step(SlabState(wet_profile),ArealSlabParameters(90000.,20.,p.wet_area_m2),30,.3,400.,100.,200*p.wet_area_m2,0.)
                wet_profile=r['state'].temperatures_C
                source+=r['water_to_debris_heat_J']+r['bare_ice_heat_J']
                ice+=r['total_ice_heat_J']
            else:
                r=water_only(ThermalState(0),ThermalParameters(0,0,100,0),30,.3,400.)
                source+=r.bare_ice_heat_J;ice+=r.bare_ice_heat_J
            if p.dry_area_m2>0:
                r=slab_step(SlabState(dry_profile),ArealSlabParameters(90000.,20.,p.dry_area_m2),30,0.,0.,0.,0.,50*p.dry_area_m2)
                dry_profile=r['state'].temperatures_C
                source+=r['atmospheric_heat_J'];ice+=r['covered_ice_heat_J']
            p=P(p.capacity_J_m2_K,p.wet_area_m2,p.dry_area_m2,wet_profile,dry_profile,p.area_reference_m2).validate()
            residual=initial+source+imported-p.stored_energy_J-ice
            self.assertLess(abs(residual)/max(1,initial+source+imported),1e-12)


if __name__=='__main__':unittest.main()
