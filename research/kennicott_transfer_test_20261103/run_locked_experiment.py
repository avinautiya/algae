"""Nominal predeclared Kennicott experiment, with separate published role gates."""
import argparse,hashlib,json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path[:0]=[str(ROOT/'kennicott_model'),str(ROOT/'algorithm'),str(ROOT/'kennicott_protocol')]
from experiment import MODEL_IDS,PRIMARY,THERMAL,INDICES,UnsupportedModel,calibrate_capacity,forward,index_doses
from outcome_loader import load_after_freeze

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):
 with Path(p).open('x')as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n')
def check(directory):
 spec=load(directory/'FROZEN_SPEC.json')
 for name,h in spec['source_hashes'].items():
  if sha(ROOT/name)!=h:raise ValueError('Frozen source changed: '+name)
 return spec

def calibrate(directory,receipt):
 spec=check(directory);truth=load_after_freeze(receipt,expected_role='train');write(directory/'TRAIN_RELEASED.json',truth)
 # Only one declared2020 target is returned by the role-gated loader.
 rows=truth['rows'] if 'rows'in truth else truth.get('released_rows')
 if rows is None:raise ValueError('Explicit source-row schema expected; never guess target field')
 if len(rows)!=1:raise ValueError('Exactly one training row required')
 row=rows[0]
 if 'observed_net_ablation_mmwe'not in row:raise ValueError('Explicit net-loss field expected; no unit guessing')
 y=row['observed_net_ablation_mmwe']
 if not math.isfinite(y)or y<=0:raise UnsupportedModel('Physical primary requires positive calibration ablation')
 plan=load(ROOT/'kennicott_model/FEATURE_ONLY_PLAN_FINAL.json');basis=load(ROOT/'kennicott_model/NOMINAL_REFERENCE_EXPOSURES_FINAL.json');models={};fits={}
 for method in [PRIMARY,'physical_homogeneous_local','physical_weather_boundary_layered']:
  try:fits[method]=dict(status='calibrated',**calibrate_capacity(plan,method,y))
  except UnsupportedModel as e:fits[method]=dict(status='refused',reason=str(e))
 fits['physical_fixed2020_layered']=fits[PRIMARY]
 for method in MODEL_IDS:
  if method in THERMAL:
   fit=fits[method]
   if fit['status']!='calibrated':models[method]=dict(status='refused',reason=fit['reason']);continue
   values=[]
   for site in ['P11','P12','P13']:
    result=forward(plan,site,method,fit['capacity_J_m3K'])
    values.append(dict(site=site,start='2023-06-08',end='2023-09-02',predicted_mmwe=result['physical_melt_mmwe'],relative_energy_residual=result['relative_energy_residual']))
   models[method]=dict(status='predicted',rows=values)
  else:
   key='homogeneous_empirical_reference_dose'if method=='homogeneous_empirical_scale'else None
   d0=basis['sites']['P_AWS_1'][key]if key else basis['sites']['P_AWS_1']['index_doses'][method]
   if d0<=0:raise ValueError('Positive training exposure required')
   values=[dict(site=s,start='2023-06-08',end='2023-09-02',predicted_mmwe=y*(basis['sites'][s][key]if key else basis['sites'][s]['index_doses'][method])/d0)for s in ['P11','P12','P13']]
   models[method]=dict(status='predicted',rows=values,global_training_gain=y/d0,physical_energy_conservation_claim=False)
 result=dict(primary=PRIMARY,models=models,calibrations=fits,training_net_ablation_loss_mmwe=y,unit='mmwe',scope='Locked retrospective conditional whole-season net-loss inference; no evaluation outcomes read',evaluation_values_opened=False)
 write(directory/'PREDICTIONS.json',result);write(directory/'PREDICTION_SEAL.json',dict(predictions_sha256=sha(directory/'PREDICTIONS.json'),frozen_spec_sha256=sha(directory/'FROZEN_SPEC.json')))
 return result
def score(directory,receipt):
 check(directory);pred=load(directory/'PREDICTIONS.json');seal=load(directory/'PREDICTION_SEAL.json')
 if sha(directory/'PREDICTIONS.json')!=seal['predictions_sha256']:raise ValueError('Prediction seal changed')
 if pred['models'][PRIMARY]['status']!='predicted':raise UnsupportedModel('Primary unscorable: preserve test outcomes')
 truth=load_after_freeze(receipt,expected_role='evaluation');rows=truth['rows']
 if [r['site']for r in rows]!=['P11','P12','P13']:raise ValueError('Frozen stake order changed')
 metrics={};refusals={}
 for name,model in pred['models'].items():
  if model['status']!='predicted':refusals[name]=model;continue
  if [r['site']for r in model['rows']]!=['P11','P12','P13']:raise ValueError('Prediction stake set mismatch')
  errors=[p['predicted_mmwe']-o['observed_net_ablation_mmwe']for p,o in zip(model['rows'],rows)]
  metrics[name]=dict(MAE_mmwe=sum(abs(x)for x in errors)/3,signed_errors_mmwe=errors,absolute_errors_mmwe=[abs(x)for x in errors],predictions=model['rows'])
 sys.path.insert(0,str(ROOT/'review'));from decision_rule import evaluate
 result=dict(primary=PRIMARY,metrics=metrics,refusals=refusals,decision=evaluate(metrics),observations=rows,predictions_sha256=sha(directory/'PREDICTIONS.json'),evaluation_receipt_sha256=sha(receipt),independent_glaciers=1,evaluation_years=1,evaluation_stakes=3,statistical_significance_established=False,scope='Locked retrospective conditional whole-season net-loss comparison; no universal or publication outcome guarantee')
 write(directory/'EVALUATION_RELEASED.json',truth);write(directory/'RESULTS.json',result);return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['calibrate','score']);p.add_argument('directory',type=Path);p.add_argument('--receipt',required=True);a=p.parse_args()
 r=calibrate(a.directory,a.receipt)if a.mode=='calibrate'else score(a.directory,a.receipt)
 print(json.dumps({'methods':len(r['models']),'primary_status':r['models'][PRIMARY]['status']}if a.mode=='calibrate'else r['decision'],indent=2))
