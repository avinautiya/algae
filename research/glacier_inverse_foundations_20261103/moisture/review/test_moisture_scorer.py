import copy,hashlib,json,tempfile,unittest
from pathlib import Path
import moisture_scorer as s

class LockedMoistureScorerTests(unittest.TestCase):
 def fixture(self):
  op=dict(quantity='native_ice_ablation_mmwe',observation_method='native_camera_917_conversion',ice_density_kg_m3=917,clock_basis='source_local_naive_no_utc_claim')
  p=dict(source_id='synthetic_Tower2',target_operator=op,heldout_targets_seen_before_prediction_freeze=False,primary_candidate=s.PRIMARY,baseline_ids=s.BASELINES,primary_metric='mean_absolute_period_ablation_error_mmwe',tie_tolerance_mmwe=1e-9,development_fit_sha256=s.EXPECTED_DEVELOPMENT_FIT_SHA256,development_parameters=s.FIXED_PARAMS.copy(),development_parameter_source='Tower1_released_periods_only',models={m:[dict(start=a,end=b,predicted_mmwe=10 if m==s.PRIMARY else 12)for a,b in zip(s.ENDPOINTS[:-1],s.ENDPOINTS[1:])]for m in s.EXPECTED})
  truth=dict(source_id=p['source_id'],target_operator=op,intervals=[dict(start=a,end=b,observed_mmwe=10)for a,b in zip(s.ENDPOINTS[:-1],s.ENDPOINTS[1:])])
  return p,truth
 def score(self,p,truth,alter=None,loader=None):
  with tempfile.TemporaryDirectory() as d:
   pred=Path(d)/'p.json';seal=Path(d)/'s.json';pred.write_text(json.dumps(p));lock=hashlib.sha256(pred.read_bytes()).hexdigest();s.write_prediction_seal(pred,seal)
   if alter:alter(pred,seal)
   return s.score_after_freeze(pred,seal,expected_prediction_sha256=lock,freeze_commit='a'*40,truth_loader=loader or (lambda:truth))
 def test_scores_entire_nine_model_family(self):
  p,t=self.fixture();r=self.score(p,t)
  self.assertEqual(len(r['model_metrics']),9);self.assertEqual(r['evaluation_calendar_periods'],4)
  self.assertTrue(r['candidate_nominally_beats_every_predeclared_baseline']);self.assertFalse(r['statistical_superiority_established'])
 def test_no_truth_opened_after_prediction_tamper(self):
  p,t=self.fixture();calls=[]
  with self.assertRaises(ValueError):self.score(p,t,lambda pred,seal:pred.write_text(pred.read_text()+' '),lambda:calls.append(1))
  self.assertEqual(calls,[])
 def test_rejects_parameter_retune_and_missing_baseline(self):
  p,_=self.fixture();p['development_parameters']['anomaly_exp_air']=5.5
  with self.assertRaises(ValueError):s.validate_predictions(p)
  p,_=self.fixture();del p['models']['calendar_exp_air']
  with self.assertRaises(ValueError):s.validate_predictions(p)
 def test_rejects_changed_camera_support(self):
  p,t=self.fixture();t['intervals'][0]['end']='2022-12-30T00:00:00'
  with self.assertRaises(ValueError):self.score(p,t)
 def test_rejects_nan_truth_no_silent_interval_drop(self):
  p,t=self.fixture();t['intervals'][1]['observed_mmwe']=float('nan')
  with self.assertRaises(ValueError):self.score(p,t)
 def test_mae_is_not_cancelling_total(self):
  p,t=self.fixture()
  for row,value in zip(p['models'][s.PRIMARY],[5,15,5,15]):row['predicted_mmwe']=value
  r=self.score(p,t);m=r['model_metrics'][s.PRIMARY]
  self.assertEqual(m['signed_total_bias_mmwe'],0);self.assertEqual(m['mean_absolute_period_ablation_error_mmwe'],5)
  self.assertFalse(r['candidate_nominally_beats_every_predeclared_baseline'])
 def test_rejects_missing_external_lock_before_truth(self):
  p,t=self.fixture()
  with tempfile.TemporaryDirectory() as d:
   a=Path(d)/'p';b=Path(d)/'s';a.write_text(json.dumps(p));s.write_prediction_seal(a,b);calls=[]
   with self.assertRaises(ValueError):s.score_after_freeze(a,b,expected_prediction_sha256='x',freeze_commit='bad',truth_loader=lambda:calls.append(1))
   self.assertEqual(calls,[])
 def test_truth_operator_density_must_match(self):
  p,t=self.fixture();t=copy.deepcopy(t);t['target_operator']['ice_density_kg_m3']=900
  with self.assertRaises(ValueError):self.score(p,t)

if __name__=='__main__':unittest.main()
