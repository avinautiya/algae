import datetime as dt,unittest
from arolla_baseline import timestamp,map_window,LOCAL
def fake():return [dict(label_local=dt.datetime(2010,8,28,h,tzinfo=LOCAL).isoformat(),model_melt_mwe_per_step=1.) for h in range(5)]
class Tests(unittest.TestCase):
 def test_explicit_local_clock(self):self.assertEqual(timestamp('04-Jul-2010 12:00:00').utcoffset(),dt.timedelta(hours=2))
 def test_label_interpretations_not_silent(self):
  a=dt.datetime(2010,8,28,1,tzinfo=LOCAL);b=a+dt.timedelta(hours=2)
  self.assertEqual(map_window(fake(),a,b,label_convention='interval_ending')['row_indices_zero_based'],[2,3]);self.assertEqual(map_window(fake(),a,b,label_convention='interval_starting')['row_indices_zero_based'],[1,2])
 def test_missing_convention_refuses(self):
  with self.assertRaises(ValueError):map_window(fake(),dt.datetime(2010,8,28,1,tzinfo=LOCAL),dt.datetime(2010,8,28,2,tzinfo=LOCAL),label_convention=None)
 def test_partial_or_naive_refuses(self):
  for a in [dt.datetime(2010,8,28,1),dt.datetime(2010,8,28,1,30,tzinfo=LOCAL)]:
   with self.assertRaises(ValueError):map_window(fake(),a,dt.datetime(2010,8,28,3,tzinfo=LOCAL),label_convention='interval_ending')
if __name__=='__main__':unittest.main()
