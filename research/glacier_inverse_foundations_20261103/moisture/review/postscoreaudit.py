"""Independent issued-prediction audit. Run ONLY after parent releases EVAL.

Never loads camera MAT files. Reads the already-authorized RESULTS artifact.
Does not import the production dose/model implementation.
"""
import argparse,bisect,csv,datetime as dt,hashlib,json,math
from pathlib import Path
from moisture_scorer import EXPECTED,FIXED_PARAMS,PRIMARY,BASELINES

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def rows(p):
 with Path(p).open()as f:return list(csv.DictReader(f))
def positive_linear_integral(u,v,duration):
 if u>=0 and v>=0:return (u+v)*duration/2
 if u<=0 and v<=0:return 0.
 return max(u,v)**2*duration/(2*abs(v-u))
def native_air_area(a,b,air):
 ts=[x[0]for x in air];i=bisect.bisect_right(ts,a)-1
 if i<0 or b>ts[-1]:raise ValueError('air extrapolation')
 total=0.
 while i+1<len(air)and air[i][0]<b:
  x,u=air[i];y,v=air[i+1];left=max(x,a);right=min(y,b)
  if right>left:
   gap=(y-x).total_seconds()
   if gap>64800:raise ValueError('air gap exceeds18h')
   vl=u+(v-u)*(left-x).total_seconds()/gap;vr=u+(v-u)*(right-x).total_seconds()/gap
   total+=positive_linear_integral(vl,vr,(right-left).total_seconds())
  i+=1
 return total

def reconstruct(features,air,temperatures,epochs,parameters,offset=0.):
 """Independent native rows + exact positive-linear air integration."""
 labels=[dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL'])for r in features]
 if any((b-a).total_seconds()!=1800 for a,b in zip(labels,labels[1:])):raise ValueError('halfhour support')
 endpoints=[labels.index(dt.datetime.fromisoformat(x))for x in epochs]
 theta=[float(r['SM3_raw_m3_m3'])+offset for r in features[:-1]]
 ref=math.fsum(theta[endpoints[0]:endpoints[1]])/(endpoints[1]-endpoints[0]);per={m:[]for m in EXPECTED}
 for j,(a,b)in enumerate(zip(labels[:-1],labels[1:])):
  ta=native_air_area(a,b,air);tr=temperatures[a.isoformat()];z=theta[j]-ref
  gradient=max((float(tr['DT4'])-float(tr['DT5']))/.11,0)*1800
  day=((a-labels[0]).total_seconds()+900)/86400-7
  q={'air_pdd':ta,'constant_rate':1800.,'buried_pdd':max(float(tr['DT1']),0)*1800,'gradient_index':gradient,
     'anomaly_exp_air':ta*math.exp(parameters['anomaly_exp_air']*z),
     'shifted_ratio_air':ta*((theta[j]+.031)/(ref+.031))**parameters['shifted_ratio_air'],
     'centered_softplus_air':ta*math.log1p(math.exp(parameters['centered_softplus_air']*z)),
     'calendar_exp_air':ta*math.exp(parameters['calendar_exp_air']*day),
     'moisture_only_exp':1800*math.exp(parameters['moisture_only_exp']*z)}
  for m in EXPECTED:per[m].append(q[m])
 return {m:[math.fsum(v[a:b])for a,b in zip(endpoints[:-1],endpoints[1:])]for m,v in per.items()}

def recalibrated_pairwise_lipschitz_per_mmice(candidate_ratios,baseline_ratios):
 """Conservative all-shared-readings bound with local TRAIN recalibration.

 Six camera endpoint errors perturb a native TRAIN total and four eval totals.
 Ratios are frozen forcing/model parameters, not refitted on camera changes.
 """
 if len(candidate_ratios)!=4 or len(baseline_ratios)!=4 or any(not math.isfinite(x)or x<0 for x in candidate_ratios+baseline_ratios):
  raise ValueError('four finite nonnegative fixed ratios per method required')
 return .917*sum(2*c+2*b+4 for c,b in zip(candidate_ratios,baseline_ratios))/4

def audit(directory,root,oldroot,expected_prediction_sha256,expected_publication_commit):
 directory=Path(directory);root=Path(root);oldroot=Path(oldroot)
 if not(directory/'RESULTS.json').exists():raise ValueError('EVAL not yet released/scored; do not access targets')
 f=read(directory/'FREEZE.json');norm=read(directory/'NORMALIZED_PREDICTIONS.json');p=read(directory/'PREDICTIONS.json');seal=read(directory/'PREDICTION_SEAL.json');r=read(directory/'RESULTS.json')
 checks=[]
 def check(name,ok):
  checks.append({'check':name,'passed':bool(ok)})
  if not ok:raise AssertionError(name)
 def close(name,a,b):check(name,math.isclose(a,b,rel_tol=2e-11,abs_tol=2e-8))
 check('published prediction hash',digest(directory/'PREDICTIONS.json')==expected_prediction_sha256==seal['prediction_sha256']==r['prediction_sha256'])
 check('published supplied commit',r['freeze_commit_supplied']==expected_publication_commit)
 check('all nine family',set(p['models'])==set(r['model_metrics'])==EXPECTED)
 check('primary unchanged',p['primary_candidate']==r['primary_candidate']==PRIMARY)
 check('development parameters unchanged',p['development_parameters']==FIXED_PARAMS)
 check('normalized artifact locked',p['normalized_predictions_sha256']==digest(directory/'NORMALIZED_PREDICTIONS.json'))
 check('physical conservation not asserted',p['physical_energy_conservation_claim']is False)
 for name,h in {**f['source_hashes'],**f['input_hashes']}.items():check('frozen '+name,digest(root/name)==h)
 for name,h in f['external_source_hashes'].items():check('external '+name,digest(oldroot/name)==h)
 feat=rows(root/'data/moisture/Tower2_CANONICAL_FEATURES.csv');rawsm={x['SOURCE_LOCAL_LABEL']:x for x in rows(oldroot/'data/decoded/sm_Tower2.csv')}
 for x in feat:close('native SM3 source row '+x['SOURCE_LOCAL_LABEL'],float(x['SM3_raw_m3_m3']),float(rawsm[x['SOURCE_LOCAL_LABEL']]['SM3']))
 air=[(dt.datetime.fromisoformat(x['SOURCE_LOCAL_LABEL']),float(x['TA_3']))for x in rows(root/'data/moisture/Tower2_NATIVE_AIR_TA3.csv')]
 temp={x['SOURCE_LOCAL_LABEL']:x for x in rows(oldroot/'data/decoded/debtemp_Tower2.csv')}
 d=reconstruct(feat,air,temp,norm['epochs'],f['parameters']);shifted=reconstruct(feat,air,temp,norm['epochs'],f['parameters'],offset=.2)
 y=p['local_train_target']['observed_ice_ablation_mmwe']
 for m in EXPECTED:
  for i,(a,b)in enumerate(zip(d[m],norm['doses'][m])):close(f'{m} independent native dose{i}',a,b)
  close(m+' local TRAIN gain',p['local_train_gains'][m],y/d[m][0])
  for j,row in enumerate(p['models'][m]):close(f'{m} prediction{j}',row['predicted_mmwe'],y*d[m][j+1]/d[m][0])
  observed=[x['observed_mmwe']for x in r['model_metrics'][m]['period_errors']]
  errors=[x['predicted_mmwe']-o for x,o in zip(p['models'][m],observed)]
  close(m+' MAE',sum(abs(x)for x in errors)/4,r['model_metrics'][m]['mean_absolute_period_ablation_error_mmwe'])
  close(m+' total bias',sum(errors),r['model_metrics'][m]['signed_total_bias_mmwe'])
 for m in ['anomaly_exp_air','moisture_only_exp','centered_softplus_air']:
  for j in range(1,5):close(f'{m} exact constant-offset calibrated ratio{j}',d[m][j]/d[m][0],shifted[m][j]/shifted[m][0])
 c=r['model_metrics'][PRIMARY]['mean_absolute_period_ablation_error_mmwe'];pairs={}
 for m in BASELINES:
  difference=c-r['model_metrics'][m]['mean_absolute_period_ablation_error_mmwe'];pair=r['primary_comparisons'][m]
  close(m+' paired MAE difference',difference,pair['candidate_minus_baseline_mean_absolute_period_error_mmwe'])
  tol=p['tie_tolerance_mmwe'];status='tie_within_frozen_tolerance'if abs(difference)<=tol else'candidate_nominal_win'if difference<0 else'baseline_nominal_win'
  check(m+' status',status==pair['status'])
  recalibration_lip=recalibrated_pairwise_lipschitz_per_mmice([x/d[PRIMARY][0]for x in d[PRIMARY][1:]],[x/d[m][0]for x in d[m][1:]])
  pairs[m]={'conditional_shared_TRAIN_and_EVAL_endpoint_recalibrated_threshold_mm_ice':abs(difference)/recalibration_lip,'recalibrated_threshold_is_sufficient_conservative_not_exact':True,'candidate_minus_baseline_MAE_mmwe':difference,'status':status,'conditional_shared_camera_endpoint_error_threshold_mm_ice':abs(difference)/(4*.917),'threshold_is_sensitivity_not_measured_precision':True}
 return {'checks':checks,'passed_checks':len(checks),'all_passed':True,'model_count':9,'primary':PRIMARY,'pairwise':pairs,
  'constant_offset_property_is_algebra_not_field_accuracy':True,'measurement_error_condition':'For FIXED issued forecasts and endpoint ice-reading errors |e|<=delta, each period target error<=2*.917*delta_mmice; pairwise MAE perturbation<=4*.917*delta_mmice. Shared endpoints do not create independence. TRAIN recalibration uncertainty is outside this fixed-forecast bound.',
  'statistical_confidence_interval':False,'historical_originality_proven':False,'universal_superiority_proven':False}

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory');p.add_argument('--root',required=True);p.add_argument('--oldroot',required=True);p.add_argument('--prediction-sha',required=True);p.add_argument('--published-commit',required=True);p.add_argument('--output',required=True);a=p.parse_args()
 result=audit(a.directory,a.root,a.oldroot,a.prediction_sha,a.published_commit)
 with Path(a.output).open('x')as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
 print(json.dumps({'all_passed':result['all_passed'],'checks':result['passed_checks'],'pairwise':result['pairwise']},indent=2))
