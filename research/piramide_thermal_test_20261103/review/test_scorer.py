import copy,json,tempfile,unittest
from pathlib import Path
from scorer import score_after_freeze,write_prediction_seal,validate_predictions

def fixture():
 rows=lambda vals:[dict(start='2022-12-16T03:00:00',end='2022-12-17T03:00:00',predicted_mmwe=vals[0]),dict(start='2022-12-17T03:00:00',end='2022-12-18T03:00:00',predicted_mmwe=vals[1])]
 p=dict(source_id='declared_synthetic_test_not_a_glacier',heldout_targets_seen_before_prediction_freeze=False,target_operator=dict(quantity='native_ice_ablation_mmwe',observation_method='synthetic camera difference',ice_density_kg_m3=917.,clock_basis='source_local_naive_no_utc_claim'),primary_candidate='candidate',baseline_ids=['rate','tie'],primary_metric='absolute_total_ablation_error_mmwe',tie_tolerance_mmwe=0.,models={'candidate':rows([20,0]),'rate':rows([6,6]),'tie':rows([8,12])})
 truth=dict(source_id=p['source_id'],target_operator=copy.deepcopy(p['target_operator']),intervals=[{k:v for k,v in r.items() if k!='predicted_mmwe'}|dict(observed_mmwe=10.) for r in rows([0,0])],uncertainty_semantics='unknown')
 return p,truth
class TestScorer(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name);self.p,self.truth=fixture();self.path=self.root/'pred.json';self.sealpath=self.root/'seal.json'
 def seal(self):
  self.path.write_text(json.dumps(self.p));return write_prediction_seal(self.path,self.sealpath)
 def score(self,seal,loader=None):return score_after_freeze(self.path,self.sealpath,expected_prediction_sha256=seal['prediction_sha256'],freeze_commit='a'*40,truth_loader=loader or (lambda:self.truth))
 def test_primary_secondary_metrics_and_ties(self):
  result=self.score(self.seal());m=result['model_metrics']['candidate']
  self.assertEqual(m['absolute_total_ablation_error_mmwe'],0);self.assertEqual(m['secondary_interval_mae_mmwe'],10)
  self.assertEqual(result['primary_comparisons']['rate']['status'],'candidate_nominal_win')
  self.assertEqual(result['primary_comparisons']['tie']['status'],'tie_within_frozen_tolerance')
  self.assertFalse(result['candidate_nominally_beats_every_predeclared_baseline']);self.assertFalse(result['statistical_superiority_established'])
 def test_bad_external_hash_before_truth_load(self):
  seal=self.seal();calls=[]
  with self.assertRaises(ValueError):score_after_freeze(self.path,self.sealpath,expected_prediction_sha256='b'*64,freeze_commit='a'*40,truth_loader=lambda:calls.append(True))
  self.assertEqual(calls,[])
 def test_postseal_prediction_change_before_truth_load(self):
  seal=self.seal();self.path.write_text(self.path.read_text()+' ');calls=[]
  with self.assertRaises(ValueError):self.score(seal,lambda:calls.append(True))
  self.assertEqual(calls,[])
 def test_scorer_change_before_truth_load(self):
  seal=self.seal();seal['scorer_sha256']='c'*64;self.sealpath.write_text(json.dumps(seal));calls=[]
  with self.assertRaises(ValueError):self.score(seal,lambda:calls.append(True))
  self.assertEqual(calls,[])
 def test_empty_baselines_or_missing_model_refused(self):
  for mutate in [lambda p:p.update(baseline_ids=[]),lambda p:p['models'].pop('tie'),lambda p:p.update(heldout_targets_seen_before_prediction_freeze=True)]:
   p=copy.deepcopy(self.p);mutate(p)
   with self.assertRaises(ValueError):validate_predictions(p)
 def test_operator_density_and_support_mismatch_refused(self):
  seal=self.seal()
  for mutate in [lambda t:t['target_operator'].update(ice_density_kg_m3=900),lambda t:t['intervals'][0].update(start='2022-12-16T04:00:00')]:
   t=copy.deepcopy(self.truth);mutate(t)
   with self.assertRaises(ValueError):self.score(seal,lambda:t)
 def test_retains_secondary_geometry_models(self):
  self.p['secondary_model_ids']=['coarse_geometry'];self.p['models']['coarse_geometry']=copy.deepcopy(self.p['models']['rate'])
  r=self.score(self.seal());self.assertEqual(len(r['model_metrics']),4);self.assertEqual(len(r['primary_comparisons']),2)
 def test_negative_camera_interval_retained_not_clipped(self):
  self.truth['intervals'][0]['observed_mmwe']=-1
  r=self.score(self.seal());self.assertEqual(r['observed_negative_intervals'],1);self.assertEqual(r['model_metrics']['rate']['observed_total_mmwe'],9)
 def test_invented_timezone_and_gaps_refused(self):
  for value in ['2022-12-16T03:00:00+00:00','2022-12-16T03:30:00']:
   p=copy.deepcopy(self.p);p['models']['rate'][0]['start']=value
   with self.assertRaises(ValueError):validate_predictions(p)
if __name__=='__main__':unittest.main()
