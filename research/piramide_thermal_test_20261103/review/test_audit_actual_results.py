import copy,csv,datetime as dt,itertools,json,tempfile,unittest
from pathlib import Path
from audit_actual_results import audit,precision_thresholds,independent_feature_dose,sha,METHODS,GEOMETRIES,HOLDS,PRIMARY,BASELINES
from scorer import write_prediction_seal,score_after_freeze

class TestActualAudit(unittest.TestCase):
 def fixture(self,root):
  start='2022-12-02T03:00:00';split='2022-12-16T03:00:00';end='2022-12-30T03:00:00';x=100.;y=80.;latent=333700.
  n=dict(start_local=start,split_local=split,end_local=end,results=[]);models={};cal={}
  for i,(method,g,h) in enumerate(itertools.product(METHODS,GEOMETRIES,HOLDS)):
   a,b=100+i,80+i;ident=f'{METHODS[method]}__{g}__{h}'
   n['results'].append(dict(method=method,geometry={'id':g},hold=h,status='predicted',no_refreezing_melt_conversion_eligible=True,training_normalized_basal_heat_K_m=a,evaluation_normalized_basal_heat_K_m=b))
   models[ident]=[dict(start=split,end=end,predicted_mmwe=x*b/a)]
   cal[ident]=dict(training_total_mmwe=x,training_total_dose=a,gain_mmwe_per_dose=x/a,implied_effective_bulk_capacity_J_m3K=x/a*latent)
  doses={}
  for ident,file,col,hold in [('positive_buried_T_DT1','debtemp_Tower1.csv','DT1','left'),('positive_air_TA3_linear_gapfill18h','trh_Tower1.csv','TA_3','linear')]:
   path=root/file;t=dt.datetime.fromisoformat(start);last=dt.datetime.fromisoformat(end)
   with path.open('w') as f:
    writer=csv.writer(f);writer.writerow(['SOURCE_LOCAL_LABEL',col])
    while t<=last:writer.writerow([t.isoformat(),5]);t+=dt.timedelta(minutes=30)
   dose=5*14*86400
   doses[ident]={part:dict(dose=dose,source_sha256=sha(path),reconstruction=hold) for part in ['train','evaluation']}
   cal[ident]=dict(training_total_mmwe=x,training_total_dose=dose,gain_mmwe_per_dose=x/dose)
   models[ident]=[dict(start=split,end=end,predicted_mmwe=x)]
  models['constant_rate']=[dict(start=split,end=end,predicted_mmwe=x)]
  normalized=root/'normalized.json';normalized.write_text(json.dumps(n))
  spec=root/'spec.json';spec.write_text(json.dumps(dict(primary_candidate=PRIMARY,baseline_ids=list(BASELINES),primary_train=[start,split],primary_evaluation=[split,end],latent_heat_J_kg=latent,tie_tolerance_mmwe=1e-6)))
  op=dict(quantity='native_ice_ablation_mmwe',clock_basis='source_local_naive_no_utc_claim',observation_method='SYNTHETIC endpoint difference',ice_density_kg_m3=917.)
  p=dict(source_id='SYNTHETIC TEST ONLY',models=models,calibrations=cal,calibration={'observed_ice_ablation_mmwe':x},feature_doses=doses,primary_candidate=PRIMARY,baseline_ids=list(BASELINES),secondary_model_ids=sorted(set(models)-{PRIMARY,*BASELINES}),target_operator=op,primary_metric='absolute_total_ablation_error_mmwe',tie_tolerance_mmwe=1e-6,heldout_targets_seen_before_prediction_freeze=False,normalized_predictions_sha256=sha(normalized),analysis_spec_sha256=sha(spec))
  predictions=root/'predictions.json';predictions.write_text(json.dumps(p));seal=root/'seal.json';sealed=write_prediction_seal(predictions,seal)
  r=score_after_freeze(predictions,seal,expected_prediction_sha256=sealed['prediction_sha256'],freeze_commit='a'*40,truth_loader=lambda:dict(source_id=p['source_id'],target_operator=op,intervals=[dict(start=split,end=end,observed_mmwe=y)],uncertainty_semantics='unknown'))
  results=root/'results.json';results.write_text(json.dumps(r));return normalized,predictions,seal,results,spec,root
 def test_release_guard_before_any_open(self):
  with self.assertRaisesRegex(ValueError,'AFTER evaluation release'):audit(*['/definitely/not/existing']*6)
 def test_full_fifteen_model_synthetic_recomputation(self):
  with tempfile.TemporaryDirectory() as d:
   paths=self.fixture(Path(d));r=audit(*paths,evaluation_released=True)
   self.assertTrue(r['all_checks_pass'],[c for c in r['checks'] if not c['pass_']]);self.assertEqual(r['models'],15)
 def test_wrong_result_metric_detected(self):
  with tempfile.TemporaryDirectory() as d:
   paths=self.fixture(Path(d));p=paths[3];r=json.loads(p.read_text());r['model_metrics'][PRIMARY]['absolute_total_ablation_error_mmwe']=123;p.write_text(json.dumps(r))
   self.assertFalse(audit(*paths,evaluation_released=True)['all_checks_pass'])
 def test_precision_threshold_matches_shared_reading_corners(self):
  x,y,fa,fb=100.,80.,.8,1.;threshold=precision_thresholds(x,y,fa,fb)['input_error_sensitivity_per_shared_reading_error_must_be_strictly_below_mmwe']
  for factor,expected in [(.99,True),(1.01,False)]:
   d=factor*threshold;wins=[]
   for e0,e1,e2 in itertools.product([-d,d],repeat=3):
    truex=x+e1-e0;truey=y+e2-e1;wins.append(abs(fa*truex-truey)<abs(fb*truex-truey))
   self.assertEqual(all(wins),expected)
 def test_piecewise_linear_zero_crossing_independent_integral(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'x.csv';p.write_text('SOURCE_LOCAL_LABEL,T\n2022-12-01T00:00:00,-2\n2022-12-01T01:00:00,2\n')
   q,_=independent_feature_dose(p,'T','2022-12-01T00:00:00','2022-12-01T01:00:00','linear',3600)
   self.assertEqual(q,1800.)
if __name__=='__main__':unittest.main()
