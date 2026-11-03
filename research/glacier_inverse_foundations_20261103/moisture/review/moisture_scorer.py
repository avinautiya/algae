"""Locked four-period Tower2 camera scoring; no fitting, selection, or frame bootstrap.
Commit identity is externally supplied, not authenticated by this utility.
"""
import datetime as dt,hashlib,itertools,json,math
from pathlib import Path
EXPECTED_DEVELOPMENT_FIT_SHA256='9043fe3e39a188eacbaf81275daba9c6931b4839a83df17827af06903deb9ebc'
PRIMARY='anomaly_exp_air'
BASELINES=['air_pdd','buried_pdd','gradient_index','constant_rate','shifted_ratio_air','centered_softplus_air','calendar_exp_air','moisture_only_exp']
EXPECTED={PRIMARY,*BASELINES}
ENDPOINTS=['2022-12-16T03:00:00','2022-12-30T03:00:00','2023-01-13T03:00:00','2023-01-27T03:00:00','2023-02-10T03:00:00']
FIXED_PARAMS={'anomaly_exp_air':5.4,'shifted_ratio_air':.34,'centered_softplus_air':7.4,'calendar_exp_air':-.0030999999999999986,'moisture_only_exp':13.300000000000004}
LOCK_KEYS=['source_id','primary_candidate','baseline_ids','target_operator','development_fit_sha256','development_parameters','development_parameter_source']

def sha(raw):return hashlib.sha256(raw).hexdigest()
def ishash(v,n=64):return isinstance(v,str) and len(v)==n and all(c in '0123456789abcdef' for c in v)
def numeric(v):return type(v) in (int,float) and math.isfinite(v)
def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False)
def support(rows,value_key,basis):
 if not isinstance(rows,list) or not rows:raise ValueError('nonempty interval block required')
 result=[];previous=None
 for r in rows:
  if not numeric(r.get(value_key)):raise ValueError('finite ablation required; no silent exclusion')
  a,b=r.get('start'),r.get('end')
  if not isinstance(a,str) or not isinstance(b,str):raise ValueError('explicit interval labels required')
  x,y=dt.datetime.fromisoformat(a),dt.datetime.fromisoformat(b)
  if basis=='source_local_naive_no_utc_claim':
   if x.tzinfo is not None or y.tzinfo is not None:raise ValueError('no invented timezone')
  elif basis=='source_documented_timezone':
   if x.tzinfo is None or y.tzinfo is None:raise ValueError('documented timezone required')
  else:raise ValueError('explicit source clock basis required')
  if not x<y or (previous is not None and a!=previous):raise ValueError('one contiguous ordered calendar block required')
  previous=b;result.append((a,b))
 return result

def validate_predictions(p):
 if p.get('heldout_targets_seen_before_prediction_freeze') is not False:raise ValueError('explicit unseen heldout declaration required')
 if not isinstance(p.get('source_id'),str) or not p['source_id'].strip():raise ValueError('source id required')
 op=p.get('target_operator',{})
 if op.get('quantity')!='native_ice_ablation_mmwe' or not op.get('observation_method'):raise ValueError('native ice camera operator required')
 if op.get('ice_density_kg_m3')!=917:raise ValueError('native source density required')
 models=p.get('models');candidate=p.get('primary_candidate');baselines=p.get('baseline_ids')
 if not isinstance(models,dict) or not models:raise ValueError('all model predictions required')
 if not isinstance(baselines,list) or not baselines or len(set(baselines))!=len(baselines) or candidate in baselines:raise ValueError('nonempty distinct baseline family required')
 secondary=p.get('secondary_model_ids',[])
 if not isinstance(secondary,list) or len(set(secondary))!=len(secondary) or set(secondary)&{candidate,*baselines}:raise ValueError('distinct prespecified secondary model ids required')
 if set(models)!={candidate,*baselines,*secondary}:raise ValueError('retain exactly the entire declared primary/secondary model family')
 if p.get('primary_metric')!='mean_absolute_period_ablation_error_mmwe':raise ValueError('fixed primary mean absolute PERIOD metric required')
 if secondary:raise ValueError('all9 models are primary or predeclared baselines; no hidden family')
 if candidate!=PRIMARY or baselines!=BASELINES or set(models)!=EXPECTED:raise ValueError('all9 frozen methods and fixed primary/baselines required')
 if not numeric(p.get('tie_tolerance_mmwe')) or p['tie_tolerance_mmwe']<0:raise ValueError('fixed tie tolerance required')
 ss=[support(r,'predicted_mmwe',op.get('clock_basis')) for r in models.values()]
 if any(s!=ss[0] for s in ss):raise ValueError('all models need identical camera support')
 if ss[0]!=list(zip(ENDPOINTS[:-1],ENDPOINTS[1:])):raise ValueError('exactly4 prespecified Tower2 periods required')
 if p.get('development_fit_sha256')!=EXPECTED_DEVELOPMENT_FIT_SHA256:raise ValueError('exact frozen development-fit hash required')
 if p.get('development_parameter_source')!='Tower1_released_periods_only':raise ValueError('Tower1-only development parameters required')
 if canonical(p.get('development_parameters'))!=canonical(FIXED_PARAMS):raise ValueError('fixed Tower1 development parameters required')
 if op.get('clock_basis')!='source_local_naive_no_utc_claim':raise ValueError('fixed source-local clock required')
 return ss[0]

def write_prediction_seal(prediction_path,seal_path):
 path=Path(seal_path)
 if path.exists():raise ValueError('refuse overwriting prediction seal')
 raw=Path(prediction_path).read_bytes();p=json.loads(raw);validate_predictions(p)
 seal=dict(version=1,prediction_sha256=sha(raw),scorer_sha256=sha(Path(__file__).read_bytes()),models=list(p['models']),secondary_model_ids=p.get('secondary_model_ids',[]))
 for key in LOCK_KEYS:seal[key]=p[key]
 path.write_text(json.dumps(seal,indent=2,allow_nan=False)+'\n');return seal

def score_after_freeze(prediction_path,seal_path,*,expected_prediction_sha256,freeze_commit,truth_loader):
 """Calls truth_loader ONLY after all local prediction/seal checks pass."""
 if not ishash(expected_prediction_sha256) or not ishash(freeze_commit,40):raise ValueError('external prediction hash and full freeze commit required')
 raw=Path(prediction_path).read_bytes();seal=json.loads(Path(seal_path).read_text())
 if sha(raw)!=expected_prediction_sha256 or seal.get('prediction_sha256')!=expected_prediction_sha256:raise ValueError('prediction differs from external lock')
 if seal.get('scorer_sha256')!=sha(Path(__file__).read_bytes()):raise ValueError('scorer changed after seal')
 p=json.loads(raw);ss=validate_predictions(p)
 for key in LOCK_KEYS:
  if seal.get(key)!=p[key]:raise ValueError('seal identity mismatch')
 if set(seal.get('models',[]))!=set(p['models']) or seal.get('secondary_model_ids',[])!=p.get('secondary_model_ids',[]):raise ValueError('model family differs from seal')
 truth=truth_loader()
 if truth.get('source_id')!=p['source_id'] or canonical(truth.get('target_operator'))!=canonical(p['target_operator']):raise ValueError('truth quantity/density/operator/source mismatch')
 if support(truth.get('intervals'),'observed_mmwe',p['target_operator']['clock_basis'])!=ss:raise ValueError('truth support differs from predictions')
 obs=[r['observed_mmwe'] for r in truth['intervals']];total=math.fsum(obs);metrics={}
 for name,rows in p['models'].items():
  values=[r['predicted_mmwe'] for r in rows];errors=[a-b for a,b in zip(values,obs)];bias=math.fsum(values)-total
  metrics[name]=dict(predicted_total_mmwe=math.fsum(values),observed_total_mmwe=total,absolute_total_ablation_error_mmwe=abs(bias),signed_total_bias_mmwe=bias,mean_absolute_period_ablation_error_mmwe=math.fsum(abs(e) for e in errors)/len(errors),secondary_absolute_total_bias_mmwe=abs(bias),secondary_interval_mae_mmwe=math.fsum(abs(e) for e in errors)/len(errors),period_errors=[dict(start=row['start'],end=row['end'],predicted_mmwe=value,observed_mmwe=observed,signed_error_mmwe=error,absolute_error_mmwe=abs(error)) for row,value,observed,error in zip(rows,values,obs,errors)],secondary_mean_interval_bias_mmwe=math.fsum(errors)/len(errors),negative_predicted_intervals=sum(v<0 for v in values),interval_count=len(values))
 candidate=p['primary_candidate'];err=metrics[candidate]['mean_absolute_period_ablation_error_mmwe'];comparisons={}
 for b in p['baseline_ids']:
  difference=err-metrics[b]['mean_absolute_period_ablation_error_mmwe']
  status='tie_within_frozen_tolerance' if abs(difference)<=p['tie_tolerance_mmwe'] else 'candidate_nominal_win' if difference<0 else 'baseline_nominal_win'
  comparisons[b]=dict(candidate_minus_baseline_mean_absolute_period_error_mmwe=difference,status=status)
 return dict(prediction_sha256=expected_prediction_sha256,freeze_commit_supplied=freeze_commit,commit_publication_verified_by_this_function=False,scorer_sha256=seal['scorer_sha256'],development_fit_sha256=p['development_fit_sha256'],development_parameters=p['development_parameters'],development_parameter_source=p['development_parameter_source'],source_id=p['source_id'],target_operator=p['target_operator'],start=ss[0][0],end=ss[-1][1],primary_candidate=candidate,model_metrics=metrics,primary_comparisons=comparisons,candidate_nominally_beats_every_predeclared_baseline=all(r['status']=='candidate_nominal_win' for r in comparisons.values()),observed_negative_intervals=sum(v<0 for v in obs),uncertainty_semantics=truth.get('uncertainty_semantics','unknown'),source_uncertainty_semantics_unknown=truth.get('uncertainty_semantics','unknown') not in {'per_reading_absolute_bound','per_interval_absolute_bound','source_standard_deviation','source_other_defined'},evaluation_calendar_periods=len(obs),serial_periods_share_observation_endpoints=True,independent_evaluation_calendar_blocks=None,frames_or_intervals_treated_as_independent=False,confidence_interval_computed=False,bootstrap_performed=False,statistical_superiority_established=False,professional_energy_balance_superiority_established=False,universal_superiority_established=False,primary_metric='mean_absolute_period_ablation_error_mmwe',tie_tolerance_mmwe=p['tie_tolerance_mmwe'],interpretation='Four locked serial Tower2 calendar periods: primary mean absolute period error; cumulative bias is secondary. Shared endpoints/one glacier do not supply four independent glaciers. All9 methods retained.')
