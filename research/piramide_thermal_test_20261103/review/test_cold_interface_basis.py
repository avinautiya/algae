"""Independent two-cell Schur-complement oracle; no field/camera access."""
import math,random,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parents[1]/'model'))
from cold_ice_column import simulate

def oracle(td,ti,ts,h,kd=1.,ki=2.,cd=2e6,ci=2e6,dd=.2,di=.2,L=333700.):
 gd=2*kd/dd;gi=2*ki/di;ad=cd*dd/h+2*gd;ai=ci*di/h+2*gi
 bd=cd*dd/h*td+gd*ts;bi=ci*di/h*ti
 constant=gd*bd/ad+gi*bi/ai
 slope=gd+gi-gd*gd/ad-gi*gi/ai
 interface=0. if constant>=0 else constant/slope
 d=(bd+gd*interface)/ad;i=(bi+gi*interface)/ai
 mass=constant*h/L if constant>=0 else 0.
 return d,i,interface,mass

def solve(td,ti,ts,h):
 return simulate(debris_conductivity_W_mK=1.,debris_heat_capacity_J_m3K=2e6,debris_depth_m=.2,initial_debris_temperature_C=[td],initial_ice_temperature_C=[ti],surface_temperature_C=[ts],interval_seconds=h,dt_seconds=h,ice_depth_m=.2,ice_conductivity_W_mK=2.,ice_heat_capacity_J_m3K=2e6,record_profiles=True)
class IndependentColdInterface(unittest.TestCase):
 def test_random_exact_two_cell_schur_oracle(self):
  rng=random.Random(53);branches=set()
  for _ in range(100):
   td,ti,ts,h=rng.uniform(-5,10),rng.uniform(-5,0),rng.uniform(-10,20),rng.uniform(30,3000)
   expected=oracle(td,ti,ts,h);r=solve(td,ti,ts,h);row=r['history'][0];branches.add(row['branch'])
   got=[r['final_debris_temperature_C'][0],r['final_ice_temperature_C'][0],row['interface_temperature_C'],r['melt_mmwe']]
   for a,b in zip(got,expected):self.assertAlmostEqual(a,b,places=10)
   self.assertLess(abs(r['relative_energy_residual']),1e-12)
  self.assertEqual(branches,{'melting','cold_interface'})
 def test_cold_steady_series_resistance(self):
  # Rs=.2+.1, Ts=-10, deep0: q=-100/3, Tinterface=-10/3.
  r=solve(-20/3,-5/3,-10,300)
  self.assertAlmostEqual(r['history'][0]['interface_temperature_C'],-10/3)
  self.assertAlmostEqual(r['history'][0]['debris_interface_heat_W_m2'],-100/3)
  self.assertEqual(r['melt_mmwe'],0.)
 def test_warm_steady_latent_units(self):
  r=solve(5,0,10,300)
  self.assertAlmostEqual(r['melt_mmwe'],50*300/333700)
  self.assertAlmostEqual(r['latent_heat_J_m2'],15000)
  self.assertEqual(r['complementarity_max_product_C_kg_m2'],0.)
 def test_cooling_creates_sensible_cold_content_not_refreezing(self):
  r=solve(0,0,-10,300)
  self.assertEqual(r['melt_mmwe'],0.)
  self.assertLess(r['final_sensible_storage_J_m2'],0)
  self.assertLess(r['final_ice_temperature_C'][0],0)
  self.assertLess(abs(r['relative_energy_residual']),1e-12)
if __name__=='__main__':unittest.main()
