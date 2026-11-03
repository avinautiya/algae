import copy,unittest
from coverage_audit import audit
def row(t):return dict(timestamp_utc=t,probe_temperature_C={'T1 -10cm':-2.,'T2 -5cm':-1.,'T3 -1cm':1.,'T4 -8cm':-2.,'T5 -1cm':1.},interface_flux_W_m2={'site1':3.,'site2':4.})
class Tests(unittest.TestCase):
 def test_calendar_not_machine_zone(self):
  r=audit([row('2015-11-30T18:15:00+00:00'),row('2015-12-28T18:15:00+00:00')]);self.assertEqual(r['support']['site1']['training']['source_records'],1);self.assertEqual(r['support']['site1']['evaluation']['source_records'],1)
 def test_values_do_not_select_favorable_support(self):
  rows=[row('2015-12-29T00:00:00+00:00')];a=audit(rows);rows[0]['interface_flux_W_m2']['site1']=-1e100;self.assertEqual(a,audit(rows))
 def test_missing_not_zero_filled(self):
  r=row('2015-12-29T00:00:00+00:00');r['interface_flux_W_m2']['site1']=None;d=audit([r]);self.assertEqual(d['support']['site1']['evaluation']['complete_records'],0)
 def test_naive_or_duplicate_refuses(self):
  for rows in [[row('2015-12-29T00:00:00')],[row('2015-12-29T00:00:00+00:00')]*2]:
   with self.assertRaises(ValueError):audit(rows)
if __name__=='__main__':unittest.main()
