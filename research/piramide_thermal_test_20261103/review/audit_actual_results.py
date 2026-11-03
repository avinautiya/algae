"""Independent AFTER-RELEASE audit of the one-block melt comparison.

No camera files are read. Explicit release authorization is required before any
prediction/calibration/result artifact is opened. No fitting or selection occurs.
"""
import argparse,csv,datetime as dt,hashlib,itertools,json,math
from pathlib import Path
METHODS={'training_temperature_forward_fit':'forward','training_derivative_OLS_NOT_source_NTLS':'derivative_OLS'}
GEOMETRIES=('small_gap_5mm','middle_gap_20mm','large_gap_40mm')
HOLDS=('left','right')
BASELINES=('derivative_OLS__middle_gap_20mm__left','constant_rate','positive_buried_T_DT1','positive_air_TA3_linear_gapfill18h')
PRIMARY='forward__middle_gap_20mm__left'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path):return json.loads(Path(path).read_text())
def checkclose(a,b):return math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-7)

def independent_feature_dose(path,column,start,end,hold,max_gap):
 """Integral by splitting at a sign crossing, rather than the production formula."""
 with Path(path).open() as stream:rows=list(csv.DictReader(stream))
 times=[dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL']) for r in rows];values=[float(r[column]) for r in rows]
 a,b=map(dt.datetime.fromisoformat,[start,end]);total=0.;covered=0.;largest=0.
 if not a<b or a<times[0] or b>times[-1]:raise ValueError('feature coverage missing')
 if any(t1<=t0 for t0,t1 in zip(times,times[1:])):raise ValueError('feature times not strictly increasing')
 if not all(math.isfinite(v) for v in values):raise ValueError('nonfinite feature')
 for i,(t0,t1) in enumerate(zip(times,times[1:])):
  left,right=max(a,t0),min(b,t1)
  if right<=left:continue
  duration=(t1-t0).total_seconds();width=(right-left).total_seconds()
  if duration>max_gap:raise ValueError('undeclared source gap')
  covered+=width;largest=max(largest,duration)
  if hold=='left':total+=max(0.,values[i])*width;continue
  if hold!='linear':raise ValueError('unknown reconstruction')
  x=values[i]+(values[i+1]-values[i])*(left-t0).total_seconds()/duration
  y=values[i]+(values[i+1]-values[i])*(right-t0).total_seconds()/duration
  if x*y<0:
   fraction=-x/(y-x)
   first=width*fraction;second=width-first
   total+=max(x,0.)*first/2+max(y,0.)*second/2
  else:total+=max(0.,x+y)*width/2
 if covered!=(b-a).total_seconds():raise ValueError('feature window incomplete')
 return total,largest

def precision_thresholds(x,y,fa,fb):
 """Conditional exact ranking thresholds; source precision remains unknown."""
 if fa==fb:return dict(identical_prediction_ratios=True,strict_ranking_possible=False)
 center=abs((fa+fb)*x-2*y)
 shared=min(center/(2*(fa+fb)+4),x/2,y/2)
 interval=min(center/(fa+fb+2),x,y)
 fixed=center/4
 return dict(identical_prediction_ratios=False,strict_ranking_possible=center>0,
  input_error_sensitivity_per_shared_reading_error_must_be_strictly_below_mmwe=shared,
  input_error_sensitivity_per_shared_reading_ice_error_must_be_strictly_below_m=shared/917,
  input_error_sensitivity_per_interval_error_must_be_strictly_below_mmwe=interval,
  fixed_issued_forecast_per_shared_reading_error_must_be_strictly_below_mmwe=fixed,
  source_precision_or_bound_known=False,achieved_robust_ranking_established=False,
  interpretation='Conditional bounds; not a standard deviation, confidence interval, or measured camera precision.')

def audit(normalized,predictions,seal,results,spec,decoded_root,*,evaluation_released=False):
 if evaluation_released is not True:raise ValueError('parent must explicitly authorize AFTER evaluation release')
 n,p,s,r,sp=[load(q) for q in [normalized,predictions,seal,results,spec]]
 checks=[]
 def test(ok,description):checks.append(dict(pass_=bool(ok),check=description))
 test(p['normalized_predictions_sha256']==sha(normalized),'normalized source hash')
 test(p['analysis_spec_sha256']==sha(spec),'analysis specification hash')
 test(s['prediction_sha256']==sha(predictions)==r['prediction_sha256'],'final prediction seal/result source hash')
 test(s['scorer_sha256']==r['scorer_sha256'],'scorer/result hash agreement')
 test(p['primary_candidate']==r['primary_candidate']==sp['primary_candidate']==PRIMARY,'predetermined primary unchanged')
 test(tuple(p['baseline_ids'])==tuple(sp['baseline_ids'])==BASELINES,'entire primary baseline family unchanged')
 test(p['target_operator']==r['target_operator'] and p['target_operator']['ice_density_kg_m3']==917,'identical native917 observation operator')
 test([n['start_local'],n['split_local']]==sp['primary_train'] and [n['split_local'],n['end_local']]==sp['primary_evaluation'],'exact predeclared training/evaluation dates')
 expected={f'{method}__{g}__{h}' for method,g,h in itertools.product(METHODS.values(),GEOMETRIES,HOLDS)}|{'constant_rate','positive_buried_T_DT1','positive_air_TA3_linear_gapfill18h'}
 test(set(p['models'])==set(r['model_metrics'])==expected and len(expected)==15,'all15 models retained')
 test(set(s['models'])==expected,'seal retains all15 models')
 x=p['calibration']['observed_ice_ablation_mmwe'];y=r['model_metrics'][PRIMARY]['observed_total_mmwe'];ratios={};independent_predictions={}
 test(math.isfinite(x) and x>0 and math.isfinite(y),'finite native training/evaluation observations')
 thermal_ids=[]
 for row in n['results']:
  if row['method'] not in METHODS:raise ValueError('unknown thermal method')
  ident=f"{METHODS[row['method']]}__{row['geometry']['id']}__{row['hold']}"
  thermal_ids.append(ident)
  a,b=row['training_normalized_basal_heat_K_m'],row['evaluation_normalized_basal_heat_K_m']
  ratio=b/a;ratios[ident]=ratio;independent_predictions[ident]=x*ratio
  c=p['calibrations'][ident]
  test(checkclose(c['training_total_mmwe'],x) and checkclose(c['training_total_dose'],a),'training-only gain inputs '+ident)
  test(checkclose(c['gain_mmwe_per_dose'],x/a),'gain=training/native dose '+ident)
  test(checkclose(c['implied_effective_bulk_capacity_J_m3K'],x/a*sp['latent_heat_J_kg']),'capacity dimensional conversion '+ident)
  test(row['status']=='predicted' and row['no_refreezing_melt_conversion_eligible'],'all frozen methods eligible '+ident)
 test(len(thermal_ids)==12 and len(set(thermal_ids))==12 and set(thermal_ids)==expected-{'constant_rate','positive_buried_T_DT1','positive_air_TA3_linear_gapfill18h'},'all12 normalized thermal scenarios uniquely retained')
 a,b,z=map(dt.datetime.fromisoformat,[n['start_local'],n['split_local'],n['end_local']])
 ratios['constant_rate']=(z-b).total_seconds()/(b-a).total_seconds();independent_predictions['constant_rate']=x*ratios['constant_rate']
 feature_meta={}
 for ident,file,column,hold,gap in [('positive_buried_T_DT1','debtemp_Tower1.csv','DT1','left',1800.),('positive_air_TA3_linear_gapfill18h','trh_Tower1.csv','TA_3','linear',64800.)]:
  path=Path(decoded_root)/file
  first,g0=independent_feature_dose(path,column,n['start_local'],n['split_local'],hold,gap)
  second,g1=independent_feature_dose(path,column,n['split_local'],n['end_local'],hold,gap)
  declared=p['feature_doses'][ident]
  test(all(declared[k]['source_sha256']==sha(path) for k in ['train','evaluation']),'PDD input source hashes '+ident)
  test(checkclose(first,declared['train']['dose']) and checkclose(second,declared['evaluation']['dose']),'independent feature integration '+ident)
  test(all(declared[k]['reconstruction']==hold for k in ['train','evaluation']),'declared source reconstruction '+ident)
  ratios[ident]=second/first;independent_predictions[ident]=x*ratios[ident]
  c=p['calibrations'][ident]
  test(checkclose(c['training_total_mmwe'],x) and checkclose(c['training_total_dose'],first) and checkclose(c['gain_mmwe_per_dose'],x/first),'training-only PDD gain '+ident)
  feature_meta[ident]=dict(train_max_gap_s=g0,evaluation_max_gap_s=g1,declared_max_gap_s=gap,native_complete_source_not_established=hold=='linear')
 for ident in expected:
  prediction=p['models'][ident]
  test(len(prediction)==1 and prediction[0]['start']==n['split_local'] and prediction[0]['end']==n['end_local'],'same one-block prediction support '+ident)
  actual=prediction[0]['predicted_mmwe'];test(checkclose(actual,independent_predictions[ident]),'Ypred=Ytrain*dose1/dose0 '+ident)
  metric=r['model_metrics'][ident]
  test(checkclose(metric['observed_total_mmwe'],y) and checkclose(metric['predicted_total_mmwe'],actual),'same observed/predicted totals '+ident)
  test(checkclose(metric['absolute_total_ablation_error_mmwe'],abs(actual-y)) and checkclose(metric['signed_total_bias_mmwe'],actual-y),'independent cumulative error and bias '+ident)
  test(checkclose(metric['secondary_interval_mae_mmwe'],abs(actual-y)),'one-block secondary MAE identity '+ident)
 primary_error=abs(independent_predictions[PRIMARY]-y);comparison={}
 for ident in BASELINES:
  d=primary_error-abs(independent_predictions[ident]-y)
  status='tie_within_frozen_tolerance' if abs(d)<=sp['tie_tolerance_mmwe'] else 'candidate_nominal_win' if d<0 else 'baseline_nominal_win'
  test(r['primary_comparisons'][ident]['status']==status,'independent nominal outcome '+ident)
  comparison[ident]=dict(status=status,error_difference_mmwe=d,precision_thresholds=precision_thresholds(x,y,ratios[PRIMARY],ratios[ident]))
 wins=all(c['status']=='candidate_nominal_win' for c in comparison.values())
 test(r['candidate_nominally_beats_every_predeclared_baseline']==wins,'nonvacuous all-baseline verdict')
 return dict(all_checks_pass=all(c['pass_'] for c in checks),checks=checks,source_results_sha256=sha(results),models=len(expected),primary_comparison=comparison,candidate_nominally_beats_all_primary_baselines=wins,feature_support=feature_meta,training_ablation_mmwe=x,evaluation_ablation_mmwe=y,evaluation_calendar_blocks=1,source_camera_precision_known=False,statistical_confidence_interval=False,independent_field_superiority_established=False,new_field_observations_or_models_fitted=False)

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 for key in ['normalized','predictions','seal','results','spec','decoded-root','out']:parser.add_argument('--'+key,type=Path,required=True)
 parser.add_argument('--evaluation-released',action='store_true');args=parser.parse_args()
 if args.out.exists():raise ValueError('refuse replacing a prior audit')
 result=audit(args.normalized,args.predictions,args.seal,args.results,args.spec,args.decoded_root,evaluation_released=args.evaluation_released)
 args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(all_checks_pass=result['all_checks_pass'],checks=len(result['checks']),nominal_all_baseline_win=result['candidate_nominally_beats_all_primary_baselines']),indent=2))
