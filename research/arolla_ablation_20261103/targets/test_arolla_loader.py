import unittest
from arolla_loader import source_time,scalar,load
from pathlib import Path
ROOT=Path(__file__).parent
class Tests(unittest.TestCase):
 def test_explicit_clock_24(self):self.assertEqual(source_time(dict(YEAR='2010',MONTH='7',DAY='31',HOUR='24',MINUTE='0')),'2010-08-01T00:00:00+02:00')
 def test_missing_snow_not_zero(self):self.assertIsNone(scalar('NaN'))
 def test_nonfinite_refuses(self):
  with self.assertRaises(ValueError):scalar('inf')
 def test_real_source_target_and_density(self):
  d=load(ROOT/'1_ARO_2010.zip',ROOT/'AROLLA_RECORD.json');self.assertEqual(d['source_ice_density_kg_m3'],900.);s=d['stake_targets'][0];self.assertEqual(s['stake'],'H1D');self.assertEqual([x['observed_ice_ablation_mm_water_equivalent'] for x in s['intervals']],[846.,315.]);self.assertTrue(d['snow']['silent_missing_to_zero_forbidden'])
 def test_source_has_complete_visit_date_envelope(self):
  d=load(ROOT/'1_ARO_2010.zip',ROOT/'AROLLA_RECORD.json');self.assertLess(d['forcing'][0]['source_label_local'],'2010-07-28');self.assertGreater(d['forcing'][-1]['source_label_local'],'2010-09-12')
if __name__=='__main__':unittest.main()
