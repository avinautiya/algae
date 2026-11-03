"""Training-only interval ice-ablation comparators on native stake measurements.

These are established empirical baselines, not a newly discovered melt law.
No heat-plate signal is converted to melt and no test target calibrates a model.
"""
import argparse,datetime as dt,hashlib,itertools,json,math
from pathlib import Path
if __package__:
 from .arolla_baseline import load as load_published
else:
 from arolla_baseline import load as load_published
ROOT=Path(__file__).resolve().parent
TZ=dt.timezone(dt.timedelta(hours=2))
MODELS=('constant_rate','positive_degree_day','net_shortwave_dose')
VISIT_HOURS=(0.,12.,24.)
LABEL_OFFSETS=(0.,.5,1.)
PUBLISHED_NAME='published_DEBCF_component_unresolved'

def finite(x):
 try:return type(x) in (float,int) and math.isfinite(x)
 except OverflowError:return False
def stamp(date,hour):return dt.datetime.combine(dt.date.fromisoformat(date),dt.time(),TZ)+dt.timedelta(hours=hour)
def bins(data,offset):
 if offset not in LABEL_OFFSETS:raise ValueError('undeclared hourly timestamp interpretation')
 out=[];last=None
 for row in data['forcing']:
  t=dt.datetime.fromisoformat(row['source_label_local'])
  if t.tzinfo is None or t.utcoffset()!=dt.timedelta(hours=2):raise ValueError('source-declared UTC+02 required')
  if last is not None and t-last!=dt.timedelta(hours=1):raise ValueError('exact hourly source support required')
  last=t
  end=t+dt.timedelta(hours=offset);out.append((end-dt.timedelta(hours=1),end,row['features']))
 return out
def doses(series,start,end):
 if not start<end:raise ValueError('positive interval required')
 need=(end-start).total_seconds();covered=0.;pdd=0.;sw=0.
 for a,b,x in series:
  seconds=max(0.,(min(b,end)-max(a,start)).total_seconds())
  if not seconds:continue
  if not all(finite(x.get(k)) for k in ['T','SWIN','SWOUT']):raise ValueError('missing forcing on target support; no silent filling')
  covered+=seconds;pdd+=max(x['T'],0.)*seconds/86400.;sw+=max(x['SWIN']-x['SWOUT'],0.)*seconds/1e6
 if abs(covered-need)>1e-6:raise ValueError('forcing does not fully cover visit-time scenario')
 return dict(constant_rate=need/86400.,positive_degree_day=pdd,net_shortwave_dose=sw)
def published_bins(rows,offset):
 if offset not in LABEL_OFFSETS:raise ValueError('undeclared model-output label interpretation')
 out=[];last=None
 for r in rows:
  t=dt.datetime.fromisoformat(r['label_local'])
  if last is not None and t-last!=dt.timedelta(hours=1):raise ValueError('published sequence not hourly')
  last=t;value=r['model_melt_mwe_per_step']
  if not finite(value) or value<0:raise ValueError('published native melt value invalid')
  end=t+dt.timedelta(hours=offset);out.append((end-dt.timedelta(hours=1),end,value))
 return out
def published_total(series,start,end):
 need=(end-start).total_seconds();covered=0.;total=0.
 if need<=0:raise ValueError('positive model window required')
 for a,b,mwe in series:
  seconds=max(0.,(min(b,end)-max(a,start)).total_seconds())
  if seconds:covered+=seconds;total+=mwe*1000.*seconds/3600.
 if abs(covered-need)>1e-6:raise ValueError('published output does not cover target interval')
 return total
def native_targets(stake):
 p=stake['points']
 if len(p)!=3 or len(stake['intervals'])!=2:raise ValueError('exact three-visit chronological design required')
 z=[x['source_melt_mm_water_equivalent'] for x in p]
 if not all(finite(v) for v in z):raise ValueError('native finite stake measurements required')
 train=z[1]-z[0];test=z[2]-z[1]
 if train<=0 or test<0:raise ValueError('positive calibration and nonnegative held-out ice ablation required')
 for i,value in enumerate([train,test]):
  if abs(value-stake['intervals'][i]['observed_ice_ablation_mm_water_equivalent'])>1e-8:raise ValueError('source interval/cumulative target disagreement')
 return train,test
def prediction(train_target,train_dose,test_dose):
 if not all(finite(x) for x in [train_target,train_dose,test_dose]) or train_target<=0 or train_dose<=0 or test_dose<0:raise ValueError('unidentifiable/invalid one-driver calibration')
 f=test_dose/train_dose;pred=train_target*f
 if not finite(pred):raise ValueError('prediction overflow')
 return pred,f,train_target/train_dose
def uncertainty(error,f,source_uncertainty_m,density):
 if not all(finite(v) and v>=0 for v in [f,source_uncertainty_m,density]):raise ValueError('uncertainty inputs invalid')
 delta=source_uncertainty_m*density
 reading=2*(f+1)*delta;interval=(f+1)*delta
 return dict(source_uncertainty_m=source_uncertainty_m,source_uncertainty_semantics='unspecified; NOT a standard deviation',
  if_per_reading_absolute_bound_residual_bound_mmwe=reading,if_per_interval_absolute_bound_residual_bound_mmwe=interval,
  compatible_if_per_reading_bound=abs(error)<=reading,compatible_if_per_interval_bound=abs(error)<=interval,
  confidence_interval=False,shared_middle_reading_propagated=True)
def check_data(data):
 if data.get('source_ice_density_kg_m3')!=900.:raise ValueError('source workbook density convention changed')
 if data.get('scope','').find('ice ablation')<0:raise ValueError('native ice-ablation contract required')
 stakes=data.get('stake_targets')
 if not isinstance(stakes,list) or {s.get('stake') for s in stakes}!={'H1D','G1D','F1D','E1D'}:raise ValueError('frozen cohort required')
 for s in stakes:native_targets(s)
 for offset in LABEL_OFFSETS:
  series=bins(data,offset)
  for s in stakes:
   dates=[p['date_local'] for p in s['points']]
   # Verify maximal visit-day envelope; missing boundaries must not disappear.
   doses(series,stamp(dates[0],0),stamp(dates[-1],24))
 return stakes
def summarize(rows):
 nominal=[r for r in rows if r['nominal']];summary={}
 for model in MODELS+(PUBLISHED_NAME,):
  rr=[r for r in nominal if r['model']==model];primary=next(r for r in rr if r['stake']=='H1D')
  all_s=[r for r in rows if r['model']==model and r['stake']=='H1D']
  summary[model]=dict(primary_nominal_prediction_mmwe=primary['predicted_mmwe'],primary_nominal_target_mmwe=primary['heldout_observed_mmwe'],primary_nominal_absolute_error_mmwe=primary['absolute_error_mmwe'],
   primary_timing_scenario_prediction_range_mmwe=[min(r['predicted_mmwe'] for r in all_s),max(r['predicted_mmwe'] for r in all_s)],
   nominal_stakes=len(rr),available_stake_nominal_mae_mmwe=sum(r['absolute_error_mmwe'] for r in rr)/len(rr),
   prediction_quantity_scope=('native_ice_ablation' if model!=PUBLISHED_NAME else 'published_component_unresolved'))
 return summary
def run(data_path,out,commit):
 if out.exists():raise ValueError('new output directory required')
 if len(commit)!=40 or any(x not in '0123456789abcdef' for x in commit):raise ValueError('published pre-execution commit required')
 freeze=json.loads((ROOT/'FREEZE.json').read_text())
 for name,sha in freeze['source_sha256'].items():
  if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=sha:raise ValueError('source changed after freeze')
 raw=data_path.read_bytes()
 if hashlib.sha256(raw).hexdigest()!=freeze['data_sha256']:raise ValueError('data changed after freeze')
 data=json.loads(raw);stakes=check_data(data);rows=[];published=load_published()
 for offset in LABEL_OFFSETS:
  series=bins(data,offset);professional=published_bins(published,offset)
  for stake in stakes:
   train_y,test_y=native_targets(stake);dates=[p['date_local'] for p in stake['points']]
   for h0,h1,h2 in itertools.product(VISIT_HOURS,repeat=3):
    t0,t1,t2=[stamp(date,h) for date,h in zip(dates,[h0,h1,h2])]
    a=doses(series,t0,t1);b=doses(series,t1,t2)
    for model in MODELS:
     pred,f,coef=prediction(train_y,a[model],b[model]);err=pred-test_y
     rows.append(dict(stake=stake['stake'],primary_stake=stake['stake']=='H1D',model=model,
      mean_label_end_offset_h=offset,visit_hours=[h0,h1,h2],nominal=(offset==0 and (h0,h1,h2)==(12.,12.,12.)),
      train_interval=[t0.isoformat(),t1.isoformat()],test_interval=[t1.isoformat(),t2.isoformat()],
      calibration_observed_mmwe=train_y,heldout_observed_mmwe=test_y,train_dose=a[model],test_dose=b[model],
      fitted_training_scale=coef,test_train_dose_ratio=f,predicted_mmwe=pred,signed_error_mmwe=err,absolute_error_mmwe=abs(err),
      relative_error_percent=100*err/test_y if test_y else None,
      uncertainty=uncertainty(err,f,data['source_stake_uncertainty_m'],data['source_ice_density_kg_m3'])))
    if stake['stake']=='H1D':
     pred=published_total(professional,t1,t2);err=pred-test_y
     rows.append(dict(stake='H1D',primary_stake=True,model=PUBLISHED_NAME,mean_label_end_offset_h=offset,visit_hours=[h0,h1,h2],nominal=(offset==0 and (h0,h1,h2)==(12.,12.,12.)),
      test_interval=[t1.isoformat(),t2.isoformat()],heldout_observed_mmwe=test_y,predicted_mmwe=pred,signed_error_mmwe=err,absolute_error_mmwe=abs(err),relative_error_percent=100*err/test_y,
      local_calibration_applied=False,published_property_independence_established=False,pure_ice_component_match_established=False,
      interpretation='As-published6cm series against native ICE ablation for context. Output ice/snow component definition unresolved; not a certified like-for-like professional skill comparison.',
      uncertainty=uncertainty(err,0.,data['source_stake_uncertainty_m'],data['source_ice_density_kg_m3'])))
 summary=summarize(rows)
 primary=[r for r in rows if r['stake']=='H1D'];robust=True
 for offset in LABEL_OFFSETS:
  for hours in itertools.product(VISIT_HOURS,repeat=3):
   r={x['model']:x for x in primary if x['mean_label_end_offset_h']==offset and tuple(x['visit_hours'])==hours}
   robust=robust and all(r['positive_degree_day']['absolute_error_mmwe']<r[m]['absolute_error_mmwe'] for m in ['constant_rate','net_shortwave_dose'])
 result=dict(preexecution_commit=commit,rows=rows,model_summary=summary,primary_PDD_strictly_better_than_both_baselines_every_declared_timing_scenario=robust,
  independent_measured_ice_ablation_target=True,calibration_and_evaluation_intervals_disjoint_except_shared_boundary_reading=True,
  published_professional_quantity_match_established=False,published_property_independence_established=False,claimed_professional_superiority=None,
  high_accuracy_established=False,independent_glacier_generalization_established=False,novel_model_claim=False,
  interpretation='Retrospective, chronologically held-out native sub-debris ICE ablation comparison. One new weather event; four spatially related stakes. Not blind, cliff/glacier-wide validation, or an award/publication guarantee.')
 out.mkdir(parents=True);(out/'RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps(dict(rows=len(rows),model_summary=summary,primary_robust_PDD_win=robust),indent=2))
def freeze(data_path):
 p=ROOT/'FREEZE.json'
 if p.exists():raise ValueError('freeze already exists; version amendments')
 raw=data_path.read_bytes();d=json.loads(raw);check_data(d)
 load_published()
 names=['arolla_melt.py','PROTOCOL.md','arolla_baseline.py','AROLLA_OUTPUT_MANIFEST.json']+[str(q.relative_to(ROOT)) for q in sorted((ROOT/'arolla_outputs').rglob('*.csv'))]
 source={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
 p.write_text(json.dumps(dict(source_sha256=source,data_sha256=hashlib.sha256(raw).hexdigest(),prediction_scores_not_computed=True,target_metadata_and_values_already_audited=True,primary_stake='H1D',scenarios_per_stake=81,models=list(MODELS),additional_primary_context_model=PUBLISHED_NAME,published_output_quantity_match_established=False),indent=2)+'\n');print('freeze ready; no prediction scored')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('data',type=Path);p.add_argument('--freeze',action='store_true');p.add_argument('--out',type=Path);p.add_argument('--commit');a=p.parse_args()
 if a.freeze:freeze(a.data)
 elif a.out and a.commit:run(a.data,a.out,a.commit)
 else:p.error('use --freeze or --out NEW --commit SHA')
