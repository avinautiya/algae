import unittest,datetime as d
from source_gate import support,valid
class Gate(unittest.TestCase):
 def row(self,t,**kw):
  r={'SOURCE_LOCAL_LABEL':t,'H_up_W_m2':'1','LE_up_W_m2':'2','pressure_Pa':'67000','gradient_Ri_at1m_approx':'.01','status':'conditional_neutral_proxy_only','wind_slope_sign_resolved_given_stated_accuracy':'True','T_slope_sign_resolved_given_stated_accuracy':'True'};r.update(kw);return r
 def test_full(self):
  a=d.datetime(2022,1,1);b=a+d.timedelta(minutes=30);self.assertTrue(support([self.row(a.isoformat()),self.row(b.isoformat())],a,b,'near_neutral')['complete'])
 def test_no_gap_bridge(self):
  a=d.datetime(2022,1,1);b=a+d.timedelta(hours=2);self.assertEqual(support([self.row(a.isoformat()),self.row(b.isoformat())],a,b,'nominal_profile')['supported_seconds'],0)
 def test_invalid_not_zeroed(self):self.assertFalse(valid(self.row('x',H_up_W_m2='',status='unsupported_nonpositive_wind_log_slope'),'nominal_profile'))
 def test_tiers(self):
  r=self.row('x',gradient_Ri_at1m_approx='.1');self.assertTrue(valid(r,'nominal_profile'));self.assertFalse(valid(r,'near_neutral'))
 def test_partial_overlap(self):
  a=d.datetime(2022,1,1);b=a+d.timedelta(minutes=30);self.assertEqual(support([self.row(a.isoformat()),self.row(b.isoformat())],a+d.timedelta(minutes=10),b,'measurement_resolved')['supported_seconds'],1200)
if __name__=='__main__':unittest.main()
