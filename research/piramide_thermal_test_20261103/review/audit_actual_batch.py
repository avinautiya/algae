"""After-release independent39-model,3-period arithmetic/source/gain audit.
No raw camera files are opened; observations come from the released result.
"""
import argparse,hashlib,json,math,sys
from pathlib import Path
from audit_actual_results import independent_feature_dose
from batch_scorer import EXPECTED,PRIMARY,BASELINES

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def close(a,b):return math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-7)
def audit(root,*,evaluation_released=False):
 if evaluation_released is not True:raise ValueError('parent must authorize future evaluation release')
 root=Path(root);base=root/'model/frozen_piramide_future_batch_final';prior=root/'model/frozen_piramide_v3';directory=base/'melt_comparison'
 n,p,s,r,pooled,f=[load(x) for x in [base/'NORMALIZED_HEAT_PREDICTIONS.json',directory/'PREDICTIONS.json',directory/'PREDICTION_SEAL.json',directory/'RESULTS.json',directory/'SECONDARY_POOLED_WITH_V3.json',base/'FREEZE.json']]
 oldp,oldn,oldr=[load(prior/x) for x in ['melt_comparison/PREDICTIONS.json','NORMALIZED_HEAT_PREDICTIONS.json','melt_comparison/RESULTS.json']]
 checks=[]
 def test(ok,label):checks.append(dict(pass_=bool(ok),check=label))
 test(n['freeze_sha256']==sha(base/'FREEZE.json'),'batch physical source freeze bound')
 test(p['normalized_predictions_sha256']==sha(base/'NORMALIZED_HEAT_PREDICTIONS.json'),'normalized dose artifact bound')
 test(s['prediction_sha256']==sha(directory/'PREDICTIONS.json')==r['prediction_sha256'],'issued prediction/result/seal hash identity')
 test(s['scorer_sha256']==r['scorer_sha256'],'scorer hash identity')
 test(p['source_fixed_gain_predictions_sha256']==n['source_fixed_gain_predictions_sha256']==sha(prior/'melt_comparison/PREDICTIONS.json'),'unchanged V3 gain source bound')
 test(p['prior_V3_result_sha256']==pooled['prior_result_sha256']==sha(prior/'melt_comparison/RESULTS.json'),'known first result retained by hash')
 test(p['primary_candidate']==n['primary_candidate']==r['primary_candidate']==PRIMARY and p['baseline_ids']==BASELINES,'unchanged primary/family')
 test(set(p['models'])==set(r['model_metrics'])==set(oldp['models'])==EXPECTED,'all39 models retained')
 test(r['primary_metric']==p['primary_metric']=='mean_absolute_period_ablation_error_mmwe','fixed3-period MAE primary')
 test(p['target_operator']==r['target_operator']==oldp['target_operator'],'same native917 operator/source footprint')
 test(p['fixed_gain_source_period']==oldp['calibration'] and p['refitted_parameters'] is False and n['gains_refitted'] is False and n['primary_changed'] is False,'TRAIN target/gains and primary not refitted')
 test(n['future_intervals']==f['config']['future_intervals'] and len(n['future_intervals'])==3,'all3 prespecified future periods')
 test(n['excluded_intervals']==f['config']['excluded_intervals']==p['unsupported_fourth_period'] and len(n['excluded_intervals'])==1,'fourth unsupported period excluded for every method before truth')
 test(f['future_target_values_opened'] is False and f['config']['future_target_values_opened'] is False and n['future_camera_values_opened'] is False,'freeze/prediction manifests predate future outcome opening')
 oldrows={f"{'forward' if x['method']=='training_temperature_forward_fit' else 'derivative_OLS'}__{x['geometry']['id']}__{x['hold']}__{x['capacity_id']}":x for x in oldn['results']}
 newrows={x['model_id']:x for x in n['results']}
 test(len(newrows)==len(n['results'])==36 and set(newrows)==set(oldrows),'36 unique unchanged thermal scenarios')
 recomputed={};ratios={};largest_ledger=0.;retreat=[]
 ytrain=oldp['calibration']['observed_ice_ablation_mmwe'];L=f['physical_config']['latent_heat_J_kg']
 for name,row in newrows.items():
  old=oldrows[name];c=oldp['calibrations'][name];C=row['physicalcapacity_J_m3K'];gain=c['gain_mmwe_per_dose']
  test(all(row[k]==old[k] for k in ['method','geometry','hold','capacity_id','physicalcapacity_J_m3K','alpha_m2_s','debris_conductivity_W_mK']),'identical thermal parameters '+name)
  test(row['source_version3_prefix_reproduced'] is True and row['no_refitting_or_test_boundary_state_reset'] is True,'continuous V3 prefix/state source certificate '+name)
  test(close(row['original_training_normalized_basal_heat_K_m'],old['training_normalized_basal_heat_K_m']) and close(row['prior_known_evaluation_normalized_basal_heat_K_m'],old['evaluation_normalized_basal_heat_K_m']),'unchanged original training/prior prefix integrals '+name)
  test(close(gain,ytrain/old['training_normalized_basal_heat_K_m']),'fixed original TRAIN gain '+name)
  test([(v['start'],v['end']) for v in row['future_blocks']]==[(v['start'],v['end']) for v in n['future_intervals']],'native period support '+name)
  recomputed[name]=[];ratios[name]=[]
  for block in row['future_blocks']:
   test(block['raw_melt_mmwe']>=0 and close(block['raw_melt_mmwe']*L,block['raw_latent_heat_J_m2']) and close(block['normalized_basal_heat_K_m']*C,block['raw_latent_heat_J_m2']),'raw LATENT/C normalization '+name+' '+block['start'])
   recomputed[name].append(gain*block['normalized_basal_heat_K_m']);ratios[name].append(block['normalized_basal_heat_K_m']/old['training_normalized_basal_heat_K_m'])
  ledger=row['physical_energy_ledger'];residual=ledger['final_sensible_storage_J_m2']-ledger['initial_sensible_storage_J_m2']+ledger['latent_heat_J_m2']-ledger['surface_heat_J_m2']+ledger['deep_heat_out_J_m2']
  test(close(residual,ledger['energy_residual_J_m2']) and row['complementarity_max_product']==0 and row['raw_component_ledger_does_not_apply_to_empirically_scaled_outputs'] is True,'unscaled physical ledger/complementarity '+name)
  largest_ledger=max(largest_ledger,abs(ledger['relative_energy_residual']));retreat.append(row['cumulative_raw_ice_equivalent_retreat_m'])
 from datetime import datetime
 duration=lambda a,b:(datetime.fromisoformat(b)-datetime.fromisoformat(a)).total_seconds()
 t0,t1=oldn['start_local'],oldn['split_local'];recomputed['constant_rate']=[ytrain*duration(v['start'],v['end'])/duration(t0,t1) for v in n['future_intervals']]
 for name,file,column,hold,gap in [('positive_buried_T_DT1','debtemp_Tower1.csv','DT1','left',1800.),('positive_air_TA3_linear_gapfill18h','trh_Tower1.csv','TA_3','linear',64800.)]:
  recomputed[name]=[]
  for i,block in enumerate(n['future_intervals']):
   path=root/'data/decoded'/file;dose,largest=independent_feature_dose(path,column,block['start'],block['end'],hold,gap);decl=p['feature_doses'][name][i]
   test(decl['source_sha256']==sha(path) and close(decl['dose'],dose) and decl['reconstruction']==hold,'independent source/PDD dose '+name+' '+str(i))
   recomputed[name].append(oldp['calibrations'][name]['gain_mmwe_per_dose']*dose)
 observations=[x['observed_mmwe'] for x in r['model_metrics'][PRIMARY]['period_errors']];metrics={};period_support=[(b['start'],b['end']) for b in n['future_intervals']]
 for name in EXPECTED:
  rr=r['model_metrics'][name];issued=p['models'][name];test([(v['start'],v['end']) for v in issued]==period_support,'issued identical period support '+name)
  test(all(close(x['predicted_mmwe'],v) for x,v in zip(issued,recomputed[name])),'independent fixed-gain prediction '+name)
  test([x['observed_mmwe'] for x in rr['period_errors']]==observations,'identical native observations '+name)
  errors=[pred-y for pred,y in zip(recomputed[name],observations)];mae=math.fsum(abs(e) for e in errors)/3;bias=math.fsum(errors)
  metrics[name]=dict(primary_mae_mmwe=mae,secondary_absolute_total_bias_mmwe=abs(bias),period_signed_errors_mmwe=errors)
  test(close(rr['mean_absolute_period_ablation_error_mmwe'],mae) and close(rr['signed_total_bias_mmwe'],bias) and close(rr['absolute_total_ablation_error_mmwe'],abs(bias)),'independent metrics '+name)
  expected_pooled=(oldr['model_metrics'][name]['absolute_total_ablation_error_mmwe']+3*mae)/4
  test(close(pooled['model_metrics'][name]['mean_absolute_four_period_error_mmwe'],expected_pooled),'secondary four-period context retains original result '+name)
 pairs={};tolerance=p['tie_tolerance_mmwe']
 for baseline in BASELINES:
  d=metrics[PRIMARY]['primary_mae_mmwe']-metrics[baseline]['primary_mae_mmwe'];status='tie_within_frozen_tolerance' if abs(d)<=tolerance else 'candidate_nominal_win' if d<0 else 'baseline_nominal_win'
  test(r['primary_comparisons'][baseline]['status']==status,'independent primary comparison '+baseline);pairs[baseline]=dict(difference_mmwe=d,status=status)
 sys.path.insert(0,str(root/'analysis'));from bounded_mae import issued_forecast_bound
 conditional={}
 for delta_ice in [0.,.0001,.0005,.001,.002,.005,.01]:
  delta_mmwe=delta_ice*917;row={}
  for baseline in BASELINES:
   value=issued_forecast_bound(recomputed[PRIMARY],recomputed[baseline],observations,delta_mmwe)
   nominal=pairs[baseline]['difference_mmwe'];expected=[nominal-4*delta_mmwe,nominal+4*delta_mmwe]
   test(all(close(a,b) for a,b in zip(value['certified_conservative_difference_interval'],expected)),'independent triangle bound '+baseline+' '+str(delta_ice))
   row[baseline]=value
  conditional[str(delta_ice)]=dict(assumed_per_reading_absolute_error_m_ice=delta_ice,assumption_is_not_a_measured_precision=True,comparisons=row)
 counts={baseline:0 for baseline in BASELINES};all_wins=0
 for name in newrows:
  if not name.startswith('forward__'):continue
  matched=name.replace('forward__','derivative_OLS__',1);mae=metrics[name]['primary_mae_mmwe'];family=[matched,*BASELINES[1:]]
  for canonical,b in zip(BASELINES,family):counts[canonical]+=mae<metrics[b]['primary_mae_mmwe']-tolerance
  all_wins+=all(mae<metrics[b]['primary_mae_mmwe']-tolerance for b in family)
 return dict(all_checks_pass=all(c['pass_'] for c in checks),checks=checks,model_count=39,calendar_period_count=3,source_results_sha256=sha(directory/'RESULTS.json'),observations_mmwe=observations,independent_metrics=metrics,primary_comparisons=pairs,descriptive18_forward_scenario_win_counts=counts,descriptive18_forward_scenarios_beating_all_four=all_wins,scenarios_are_not_probabilities_or_independent_samples=True,conditional_fixed_issued_forecast_bounds=conditional,source_camera_precision_known=False,conditional_bounds_are_not_confidence_intervals=True,max_raw_physical_relative_energy_residual=largest_ledger,cumulative_raw_ice_equivalent_retreat_m_range=[min(retreat),max(retreat)],full_physical_trace_not_duplicated_by_this_audit=True,raw_camera_files_opened=False,statistical_or_professional_superiority_established=False)
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--evaluation-released',action='store_true');args=parser.parse_args()
 if args.out.exists():raise ValueError('refuse overwriting audit')
 result=audit(args.root,evaluation_released=args.evaluation_released);args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ['all_checks_pass','model_count','calendar_period_count','primary_comparisons','descriptive18_forward_scenario_win_counts','descriptive18_forward_scenarios_beating_all_four']},indent=2))
