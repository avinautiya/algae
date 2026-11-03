import copy,unittest
from paired_skill import score
UNIT='mm_water_equivalent_per_day'
def rows():return [dict(key='A/1',site='A',observed=1.,predictions=dict(candidate=2.,baseline=1.),unit=UNIT),dict(key='B/1',site='B',observed=2.,predictions=dict(candidate=2.,baseline=0.),unit=UNIT)]
def run(r):return score(r,candidate='candidate',models=['candidate','baseline'],unit=UNIT,expected_keys=['A/1','B/1'])
class Tests(unittest.TestCase):
 def test_absent_support_identity_refuses(self):
  with self.assertRaises(ValueError):score([dict(key=None,site='A',observed=1.,predictions=dict(candidate=1.,baseline=2.),unit=UNIT)],candidate='candidate',models=['candidate','baseline'],unit=UNIT,expected_keys=[None])
 def test_nonobject_record_refuses(self):
  with self.assertRaises(ValueError):run(['not a record'])
 def test_no_baseline_is_not_a_win(self):
  with self.assertRaises(ValueError):score([dict(key='A',site='A',observed=1.,predictions=dict(candidate=1.),unit=UNIT)],candidate='candidate',models=['candidate'],unit=UNIT,expected_keys=['A'])
 def test_pooled_win_is_not_consistent_site_win(self):
  s=run(rows());self.assertTrue(s['candidate_lower_pooled_mae_than_every_baseline']);self.assertFalse(s['candidate_lower_mae_at_every_site_than_every_baseline'])
 def test_hand_calculation(self):
  r=run(rows());self.assertEqual(r['metrics']['candidate']['mae'],.5);self.assertEqual(r['metrics']['candidate']['bias'],.5);self.assertEqual(r['day_win_counts']['baseline']['candidate_lower_error'],1);self.assertEqual(r['day_win_counts']['baseline']['candidate_higher_error'],1)
 def test_missing_record_refuses(self):
  with self.assertRaises(ValueError):run(rows()[:1])
 def test_missing_model_refuses(self):
  r=rows();del r[0]['predictions']['baseline']
  with self.assertRaises(ValueError):run(r)
 def test_nan_and_bool_refuse(self):
  for v in [float('nan'),True,10**1000]:
   r=rows();r[0]['observed']=v
   with self.assertRaises(ValueError):run(r)
 def test_duplicate_refuses(self):
  with self.assertRaises(ValueError):run([rows()[0],rows()[0]])
 def test_units_refuse(self):
  r=rows();r[0]['unit']='metres'
  with self.assertRaises(ValueError):run(r)
 def test_signed_error_retained(self):
  r=rows();r[0]['observed']=-1.;s=run(r);self.assertEqual(s['metrics']['candidate']['bias'],1.5)
if __name__=='__main__':unittest.main()
