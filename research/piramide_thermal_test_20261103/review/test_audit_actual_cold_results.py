import copy,csv,datetime as dt,itertools,json,tempfile,unittest
from pathlib import Path
from audit_actual_cold_results import audit,precision_thresholds,independent_feature_dose,sha,METHODS,GEOMETRIES,HOLDS,PRIMARY,BASELINES
from scorer import write_prediction_seal,score_after_freeze

class TestColdActualAudit(unittest.TestCase):
 def fixture(self,root):
  start='2022-12-02T03:00:00';split='2022-12-16T03:00:00';end='2022-12-30T03:00:00';x=100.;y=80.;latent=333700.
  n=dict(start_local=start,split_local=split,end_local=end,results=[],raw_component_uses_cold_ice_complementarity=True,primary_capacity_id='proxy',camera_outcome_fields_loaded=False,calibration_melt_gain_fitted=False);models={};cal={}
  for i,(method,g,h,capacity) in enumerate(itertools.product(METHODS,GEOMETRIES,HOLDS,['dry','proxy','saturated'])):
   a,b=100+i,80+i;ident=f'{METHODS[method]}__{g}__{h}__{capacity}';C={'dry':1e6,'proxy':1.5e6,'saturated':2e6}[capacity]
   trace=[dict(start_local=start,end_local=split,normalized_basal_heat_K_m=a,raw_melt_mmwe=a*C/latent,raw_latent_heat_J_m2=a*C),dict(start_local=split,end_local=end,normalized_basal_heat_K_m=b,raw_melt_mmwe=b*C/latent,raw_latent_heat_J_m2=b*C)]
   n['results'].append(dict(method=method,geometry={'id':g},hold=h,capacity_id=capacity,physicalcapacity_J_m3K=C,capacity_J_m3K=C,status='predicted',no_refreezing_melt_conversion_eligible=True,training_normalized_basal_heat_K_m=a,evaluation_normalized_basal_heat_K_m=b,normalized_heat_trace=trace,raw_training_melt_mmwe=a*C/latent,raw_evaluation_melt_mmwe=b*C/latent,physical_energy_ledger=dict(initial_sensible_storage_J_m2=0,final_sensible_storage_J_m2=0,latent_heat_J_m2=(a+b)*C,surface_heat_J_m2=(a+b)*C,deep_heat_out_J_m2=0,energy_residual_J_m2=0),ledger_applies_to_unscaled_raw_component_only=True,empirical_output_gain_is_not_physical_heat_capacity=True,negative_latent_melt_substeps=0,interface_complementarity_max_product=0,alpha_refitted=False,wet_debris_phase_change_included=False,fixed_geometry_no_ice_thermal_advection=True))
   models[ident]=[dict(start=split,end=end,predicted_mmwe=x*b/a)]
   cal[ident]=dict(training_total_mmwe=x,training_total_dose=a,gain_mmwe_per_dose=x/a,scaled_capacity_diagnostic_NOT_physical_capacity_J_m3K=x/a*latent,postprocessing_gain_relative_to_physical_component=x/a*latent/C,gain_is_not_identified_physical_heat_capacity=True)
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
  physical=root/'physical_freeze.json';physical.write_text(json.dumps(dict(config=dict(capacities=[dict(id=c,capacity_J_m3K=C) for c,C in [('dry',1e6),('proxy',1.5e6),('saturated',2e6)]],cold_ice_initial_temperature_C=0.,ice_depth_m=2.))))
  n['freeze_sha256']=sha(physical)
  normalized=root/'normalized.json';normalized.write_text(json.dumps(n))
  spec=root/'spec.json';spec.write_text(json.dumps(dict(primary_candidate=PRIMARY,baseline_ids=list(BASELINES),primary_train=[start,split],primary_evaluation=[split,end],latent_heat_J_kg=latent,tie_tolerance_mmwe=1e-6)))
  op=dict(quantity='native_ice_ablation_mmwe',clock_basis='source_local_naive_no_utc_claim',observation_method='SYNTHETIC endpoint difference',ice_density_kg_m3=917.)
  p=dict(source_id='SYNTHETIC TEST ONLY',models=models,calibrations=cal,calibration={'observed_ice_ablation_mmwe':x},feature_doses=doses,primary_candidate=PRIMARY,baseline_ids=list(BASELINES),secondary_model_ids=sorted(set(models)-{PRIMARY,*BASELINES}),target_operator=op,primary_metric='absolute_total_ablation_error_mmwe',tie_tolerance_mmwe=1e-6,heldout_targets_seen_before_prediction_freeze=False,normalized_predictions_sha256=sha(normalized),analysis_spec_sha256=sha(spec))
  predictions=root/'predictions.json';predictions.write_text(json.dumps(p));seal=root/'seal.json';sealed=write_prediction_seal(predictions,seal)
  r=score_after_freeze(predictions,seal,expected_prediction_sha256=sealed['prediction_sha256'],freeze_commit='a'*40,truth_loader=lambda:dict(source_id=p['source_id'],target_operator=op,intervals=[dict(start=split,end=end,observed_mmwe=y)],uncertainty_semantics='unknown'))
  results=root/'results.json';results.write_text(json.dumps(r));return normalized,predictions,seal,results,spec,root
 def test_release_guard_before_any_open(self):
  with self.assertRaisesRegex(ValueError,'AFTER evaluation release'):audit(*['/definitely/not/existing']*6)
 def test_full_thirtynine_model_synthetic_recomputation(self):
  with tempfile.TemporaryDirectory() as d:
   paths=self.fixture(Path(d));r=audit(*paths,evaluation_released=True,physical_freeze=Path(d)/'physical_freeze.json')
   self.assertTrue(r['all_checks_pass'],[c for c in r['checks'] if not c['pass_']]);self.assertEqual(r['models'],39)
 def test_wrong_result_metric_detected(self):
  with tempfile.TemporaryDirectory() as d:
   paths=self.fixture(Path(d));p=paths[3];r=json.loads(p.read_text());r['model_metrics'][PRIMARY]['absolute_total_ablation_error_mmwe']=123;p.write_text(json.dumps(r))
   self.assertFalse(audit(*paths,evaluation_released=True,physical_freeze=Path(d)/'physical_freeze.json')['all_checks_pass'])
 def test_raw_latent_ledger_corruption_detected(self):
  with tempfile.TemporaryDirectory() as d:
   paths=self.fixture(Path(d));p=paths[0];n=json.loads(p.read_text());n['results'][0]['normalized_heat_trace'][0]['raw_latent_heat_J_m2']+=100;p.write_text(json.dumps(n))
   result=audit(*paths,evaluation_released=True,physical_freeze=Path(d)/'physical_freeze.json')
   self.assertFalse(result['all_checks_pass'])
   self.assertTrue(any(not c['pass_'] and 'LATENT normalization' in c['check'] for c in result['checks']))
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
