import math,unittest
from column import simulate

def scenario(**kw):
 base={'surface_temperature_C':[5.]*24,'interval_seconds':300.,'debris_depth_m':.2,'debris_alpha_cells_m2_s':[1e-6]*8,'debris_capacity_J_m3K':1.5e6,'initial_debris_temperature_C':[5*(1-(i+.5)/8)for i in range(8)],'initial_ice_temperature_C':[0.]*10,'ice_depth_m':1.}
 base.update(kw);return simulate(**base)
class ColumnTests(unittest.TestCase):
 def test_uniform_exact_steady_resistance_and_ledger(self):
  r=scenario();expected=(1.5*5/.2)*24*300/333700
  self.assertAlmostEqual(r['physical_melt_mmwe'],expected,places=11)
  self.assertLess(abs(r['relative_energy_residual']),1e-12)
 def test_heterogeneous_exact_series_resistance(self):
  k=[1.,1.,2.,2.];dz=.2/4;q=5/sum(dz/v for v in k);T=[];res=0
  for v in k:T.append(5-q*(res+dz/(2*v)));res+=dz/v
  r=scenario(debris_alpha_cells_m2_s=[v/1.5e6 for v in k],initial_debris_temperature_C=T)
  self.assertAlmostEqual(r['physical_melt_mmwe'],q*24*300/333700,places=11)
  self.assertLess(abs(r['relative_energy_residual']),1e-12)
 def test_cooling_stores_negative_heat_without_imaginary_melt(self):
  r=scenario(surface_temperature_C=[-5.]*72,initial_debris_temperature_C=[0.]*8)
  self.assertEqual(r['physical_melt_mmwe'],0.);self.assertGreater(r['cold_interface_steps'],0)
  self.assertTrue(all(x['interface_T_C']<=0 and x['melt_mmwe']==0 for x in r['trace']))
  self.assertLess(abs(r['relative_energy_residual']),1e-12)
 def test_cold_initial_ice_reduces_melt_and_closes_ledger(self):
  warm=scenario(surface_temperature_C=[5.]*288)
  cold=scenario(surface_temperature_C=[5.]*288,initial_ice_temperature_C=[-1.]*10)
  self.assertLess(cold['physical_melt_mmwe'],warm['physical_melt_mmwe'])
  self.assertLess(abs(cold['relative_energy_residual']),1e-11)
 def test_negative_material_and_ice_above_melt_refused(self):
  with self.assertRaises(ValueError):scenario(initial_ice_temperature_C=[1.]*10)
  with self.assertRaises(ValueError):scenario(debris_alpha_cells_m2_s=[-1.]*8)
if __name__=='__main__':unittest.main()
