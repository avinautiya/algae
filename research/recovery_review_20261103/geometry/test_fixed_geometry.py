import unittest
import numpy as np
from fixed_geometry import FixedGeometrySim
from vendor.swe2d_v4_pinned import _bc
class Tests(unittest.TestCase):
 def test_identical_physical_dimensions_and_areas(self):
  descriptors=[]
  for dx in [.25,.125,.0625]:
   s=FixedGeometrySim(dx=dx,min_lane_cells=1);descriptors.append(s.geometry_contract)
   self.assertAlmostEqual(s.ny*dx,.75);self.assertAlmostEqual(s.nx*dx,1.25);self.assertAlmostEqual(s.divider_cells*dx,.25)
   self.assertAlmostEqual(s.toe.sum()*dx**2,.125);self.assertAlmostEqual(s.byp.sum()*dx**2,.125);self.assertAlmostEqual(s.wall.sum()*dx**2,.1875)
   self.assertEqual(s.divider_cells,round(.25/dx));self.assertFalse(np.any(s.wall&s.ice));self.assertAlmostEqual(s.Q/(s.ny*s.dx),.07/.75)
  self.assertEqual(descriptors[0],descriptors[1]);self.assertEqual(descriptors[1],descriptors[2])
 def test_ghosts_and_initial_ledger(self):
  s=FixedGeometrySim();self.assertFalse(s.ice[0].any() or s.ice[-1].any() or s.ice[:,0].any() or s.ice[:,-1].any())
  np.testing.assert_array_equal(s.zb[:,0],s.zb[:,1]);np.testing.assert_array_equal(s.zb[:,-1],s.zb[:,-2])
  self.assertEqual(s.source_log,[]);self.assertEqual(s.hydraulic_window_log,[]);self.assertEqual(s.physical_time_s,0.)
  self.assertAlmostEqual(s.combined_ledgers()['combined_ice_water_mass_residual_kg'],0.);self.assertTrue(all(r['kind']=='initial_condition_definition' for r in s.reset_log))
 def test_reject_non_aligned_geometry(self):
  for kw in [dict(dx=.1),dict(divider_width_m=.2),dict(Ltail=.2),dict(dx=False)]:
   with self.assertRaises(ValueError):FixedGeometrySim(**kw)
 def test_route_slices_and_wall_initial_state(self):
  s=FixedGeometrySim(dx=.0625);self.assertEqual(s.byp_cols().start,9);self.assertEqual(s.byp_cols().stop,13)
  self.assertTrue(np.all(s.h[s.wall]==0));self.assertTrue(np.all(s.hT[s.wall]==0));self.assertFalse(np.any(s.toe&s.byp))
 def test_short_physical_step_no_wall_melt(self):
  s=FixedGeometrySim(dx=.0625);z=s.zb.copy();r=s.step_physical(.01)
  np.testing.assert_array_equal(s.zb[s.wall],z[s.wall]);self.assertFalse(s.ice[s.wall].any());self.assertTrue(r['ledgers']['physical_split_only'])
  self.assertLess(abs(r['ledgers']['combined_ice_water_mass_residual_kg']),1e-9);self.assertLess(abs(r['ledgers']['physical_split_thermal_latent_residual_J']),1e-5)
 def test_post_evolution_reset_exposed_not_erased(self):
  s=FixedGeometrySim();s.step_physical(.01);s.reset_state(.16);self.assertFalse(s.combined_ledgers()['physical_split_only']);self.assertNotEqual(s.combined_ledgers()['explicit_reset_water_m3'],0.)
 def test_minimum_lane_guard_is_not_convergence_certificate(self):
  with self.assertRaises(ValueError):FixedGeometrySim(dx=.125)
  s=FixedGeometrySim();self.assertEqual(s.mesh_contract['lane_cells'],4);self.assertFalse(s.mesh_contract['convergence_established'])
 def test_actual_inlet_boundary_integral_fixed_across_meshes(self):
  for dx in [.25,.125,.0625]:
   s=FixedGeometrySim(dx=dx,min_lane_cells=1)
   _bc(s.h,s.hu,s.hv,s.hT,s.nx,s.ny,s.Q/(s.ny*s.dx),s.Tw,s.closed_out)
   self.assertAlmostEqual(float(s.hu[0,1:-1].sum()*dx),s.Q,places=14)
   np.testing.assert_allclose(s.hT[0,1:-1],s.h[0,1:-1]*s.Tw)
   self.assertAlmostEqual(s.combined_ledgers()['combined_ice_water_mass_residual_kg'],0.)
if __name__=='__main__':unittest.main()
