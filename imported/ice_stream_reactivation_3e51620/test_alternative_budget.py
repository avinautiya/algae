import unittest,math
from alternative_budget import capacity,audit
class Tests(unittest.TestCase):
 def test_polygon(self):
  w,h,t,b,p=2.,3.,1.2,.6,.6
  x=-h/math.tan(t);pts=[(0.,0.),(x,h),(x+h/math.tan(b),0.)]
  area=abs(sum(pts[i][0]*pts[(i+1)%3][1]-pts[(i+1)%3][0]*pts[i][1] for i in range(3)))/2
  self.assertAlmostEqual(capacity(w,h,t,b,p),area*w*p)
 def test_bound_and_duality(self):
  d=dict(width=2.,height=3.,mantle_solid_m=.1,englacial_solid_fraction=0.,stock_capacity=1.,initial_stock=0.,external_input=0.,export_upper=0.,energy=1e20,retreat=5.)
  self.assertEqual(audit(**d)['status'],'conditionally feasible')
  d['retreat']=6.;a=audit(**d);self.assertEqual(a['status'],'incompatible')
  d['export_upper']=a['required_export_m3'];self.assertEqual(audit(**d)['status'],'conditionally feasible')
 def test_energy(self):
  a=audit(width=1.,height=1.,mantle_solid_m=0.,englacial_solid_fraction=0.,stock_capacity=1.,initial_stock=0.,external_input=0.,export_upper=0.,energy=0.,retreat=1.)
  self.assertEqual(a['status'],'incompatible')
 def test_invalid(self):
  with self.assertRaises(ValueError):capacity(1.,1.,.5,.6,.6)
if __name__=='__main__':unittest.main()
