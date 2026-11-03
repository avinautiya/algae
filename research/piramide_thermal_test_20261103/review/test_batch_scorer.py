import copy,json,math,tempfile,unittest
from pathlib import Path
from batch_scorer import EXPECTED,PRIMARY,BASELINES,validate_predictions,write_prediction_seal,score_after_freeze,pooled_with_prior_result,sha

def fixture():
 times=['2022-12-30T00:00:00','2023-01-13T03:00:00','2023-01-27T03:00:00','2023-02-10T03:00:00']
 rows=lambda vals:[dict(start=times[i],end=times[i+1],predicted_mmwe=v) for i,v in enumerate(vals)]
 models={m:rows([12,12,12]) for m in EXPECTED};models[PRIMARY]=rows([20,0,10]);models[BASELINES[-1]]=rows([0,20,10])
 p=dict(source_id='SYNTHETIC NO REAL CAMERA VALUES',heldout_targets_seen_before_prediction_freeze=False,target_operator=dict(quantity='native_ice_ablation_mmwe',observation_method='synthetic adjacent endpoint difference',ice_density_kg_m3=917.,clock_basis='source_local_naive_no_utc_claim'),primary_candidate=PRIMARY,baseline_ids=list(BASELINES),secondary_model_ids=sorted(EXPECTED-{PRIMARY,*BASELINES}),primary_metric='mean_absolute_period_ablation_error_mmwe',tie_tolerance_mmwe=1e-6,models=models)
 truth=dict(source_id=p['source_id'],target_operator=copy.deepcopy(p['target_operator']),intervals=[dict(start=times[i],end=times[i+1],observed_mmwe=10.) for i in range(3)],uncertainty_semantics='unknown')
 return p,truth
class TestBatchScorer(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name);self.p,self.truth=fixture();self.pred=self.root/'pred.json';self.seal=self.root/'seal.json'
 def publish(self):self.pred.write_text(json.dumps(self.p));return write_prediction_seal(self.pred,self.seal)
 def score(self,s,loader=None):return score_after_freeze(self.pred,self.seal,expected_prediction_sha256=s['prediction_sha256'],freeze_commit='b'*40,truth_loader=loader or (lambda:self.truth))
 def test_primary_mae_cannot_be_replaced_by_cancelled_total_bias(self):
  r=self.score(self.publish());m=r['model_metrics'][PRIMARY]
  self.assertAlmostEqual(m['mean_absolute_period_ablation_error_mmwe'],20/3);self.assertEqual(m['absolute_total_ablation_error_mmwe'],0.)
  self.assertEqual(r['primary_comparisons'][BASELINES[0]]['status'],'baseline_nominal_win')
  self.assertEqual(r['primary_comparisons'][BASELINES[-1]]['status'],'tie_within_frozen_tolerance')
  self.assertFalse(r['candidate_nominally_beats_every_predeclared_baseline']);self.assertEqual(len(r['model_metrics']),39)
  self.assertIsNone(r['independent_evaluation_calendar_blocks']);self.assertTrue(r['serial_periods_share_observation_endpoints']);self.assertFalse(r['confidence_interval_computed'])
 def test_no_empty_family_missing_scenario_or_changed_primary(self):
  for mutate in [lambda p:p.update(baseline_ids=[]),lambda p:p['models'].pop(next(iter(p['secondary_model_ids']))),lambda p:p.update(primary_candidate='constant_rate'),lambda p:p.update(primary_metric='absolute_total_ablation_error_mmwe')]:
   p=copy.deepcopy(self.p);mutate(p)
   with self.assertRaises(ValueError):validate_predictions(p)
 def test_exactly_three_periods_all_models(self):
  for model in self.p['models']:self.p['models'][model].pop()
  with self.assertRaises(ValueError):validate_predictions(self.p)
 def test_corrupt_artifact_does_not_call_future_truth_loader(self):
  s=self.publish();self.pred.write_text(self.pred.read_text()+' ');calls=[]
  with self.assertRaises(ValueError):self.score(s,lambda:calls.append(True))
  self.assertEqual(calls,[])
 def test_wrong_scorer_hash_does_not_call_loader(self):
  s=self.publish();s['scorer_sha256']='d'*64;self.seal.write_text(json.dumps(s));calls=[]
  with self.assertRaises(ValueError):self.score(s,lambda:calls.append(True))
  self.assertEqual(calls,[])
 def test_observation_operator_and_interval_match(self):
  s=self.publish();self.truth['target_operator']['ice_density_kg_m3']=900.
  with self.assertRaises(ValueError):self.score(s)
 def test_missing_value_refuses_not_zero_imputation(self):
  s=self.publish();self.truth['intervals'][1]['observed_mmwe']=float('nan')
  with self.assertRaises(ValueError):self.score(s)
 def test_pooled_four_period_context_retains_first_loss(self):
  r=self.score(self.publish());prior=dict(source_id=r['source_id'],target_operator=r['target_operator'],end=r['start'],primary_candidate=PRIMARY,model_metrics={m:dict(interval_count=1,absolute_total_ablation_error_mmwe=4.,signed_total_bias_mmwe=4.) for m in EXPECTED})
  path=self.root/'prior.json';raw=json.dumps(prior).encode();path.write_bytes(raw)
  pooled=pooled_with_prior_result(r,path,expected_prior_result_sha256=sha(raw))
  self.assertAlmostEqual(pooled['model_metrics'][PRIMARY]['mean_absolute_four_period_error_mmwe'],6.)
  self.assertFalse(pooled['replaces_three_new_period_primary_metric']);self.assertTrue(pooled['includes_prior_v3_first_period'])
  with self.assertRaises(ValueError):pooled_with_prior_result(r,path,expected_prior_result_sha256='f'*64)
if __name__=='__main__':unittest.main()
