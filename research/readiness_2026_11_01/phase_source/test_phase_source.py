import unittest
import numpy as np
from phase_source import apply,route_source_ratio

class SourceTests(unittest.TestCase):
    def state(self,**kw):
        c=dict(bed_m=np.array([0.,-1.]),depth_m=np.array([.1,.2]),momentum_x_m2_s=np.array([.03,.04]),momentum_y_m2_s=np.array([.01,0.]),depth_temperature_m_K=np.array([.1,.4]),ice_melt_depth_m=np.array([.02,.01]),cell_area_m2=np.array([2.,3.]));c.update(kw);return c
    def test_phase_mass_independent_integrals(self):
        r=apply(**self.state());self.assertAlmostEqual(r.water_volume_added_m3,.917*.07);self.assertAlmostEqual(r.ice_volume_lost_m3,.07);self.assertLess(abs(r.phase_mass_residual_kg),1e-12)
    def test_density_jump_surface(self):
        c=self.state();r=apply(**c);np.testing.assert_allclose(r.bed_m+r.depth_m,np.array([.1,-.8])-.083*np.array([.02,.01]),atol=1e-15)
    def test_equal_density_surface_fixed(self):
        c=self.state();r=apply(**c,rho_i_kg_m3=1000.);np.testing.assert_allclose(r.bed_m+r.depth_m,c['bed_m']+c['depth_m'],atol=1e-15)
    def test_zero_heat_and_momentum_birth_dilution(self):
        c=self.state();r=apply(**c);np.testing.assert_array_equal(r.depth_temperature_m_K,c['depth_temperature_m_K']);np.testing.assert_array_equal(r.momentum_x_m2_s,c['momentum_x_m2_s']);self.assertTrue(np.all(r.depth_temperature_m_K/r.depth_m<c['depth_temperature_m_K']/c['depth_m']))
    def test_dry_creation_zero_temperature(self):
        c=self.state(depth_m=np.zeros(2),momentum_x_m2_s=np.zeros(2),momentum_y_m2_s=np.zeros(2),depth_temperature_m_K=np.zeros(2));r=apply(**c);self.assertTrue(np.all(r.depth_m>0));self.assertEqual(float(r.depth_temperature_m_K.sum()),0.)
    def test_no_mutation_and_zero_melt(self):
        c=self.state(ice_melt_depth_m=np.zeros(2));copies={k:v.copy() for k,v in c.items()};r=apply(**c);r.bed_m[0]=999
        for k,v in c.items():np.testing.assert_array_equal(v,copies[k])
    def test_refusals(self):
        for kw in [dict(ice_melt_depth_m=np.array([-.1,0.])),dict(depth_m=np.array([0.,.2])),dict(cell_area_m2=0.),dict(bed_m=np.array([np.nan,0.]))]:
            with self.assertRaises(ValueError):apply(**self.state(**kw))
    def test_route_local_zero_not_global_ratio(self):
        self.assertIsNone(route_source_ratio(generated_water_m3=.001,interval_s=100.,route_discharge_m3_s=0.));self.assertAlmostEqual(route_source_ratio(generated_water_m3=.001,interval_s=100.,route_discharge_m3_s=.0001),.1)
if __name__=='__main__':unittest.main()
