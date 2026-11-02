"""Independent standard-library arithmetic/truth audit; never changes predictions."""
import csv,datetime as dt,hashlib,json,math,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parent
FIELD=ROOT.parent
MODELS=['diffusion','quasisteady','lag_filter','ridge']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def metric(truth,pred):
 e=[p-y for y,p in zip(truth,pred)];a=list(map(abs,e));n=len(e)
 return dict(n=n,MAE_C=math.fsum(a)/n,MSE_C=math.fsum(x*x for x in e)/n,RMSE_C=math.sqrt(math.fsum(x*x for x in e)/n),bias_C=math.fsum(e)/n,
 fraction_within_0_2C=sum(x<=.2 for x in a)/n,fraction_within_0_5C=sum(x<=.5 for x in a)/n,fraction_within_1C=sum(x<=1 for x in a)/n)
def run():
 result=json.loads((ROOT/'results/FIELD_TEMPERATURE_RESULTS.json').read_text());diffs=[];profiles=[];allids=set()
 for p in result['profile_results']:
  path=ROOT/'results'/p['prediction_csv'];assert sha(path)==p['prediction_sha256']
  with path.open() as f:rows=list(csv.DictReader(f))
  rawpath=FIELD/'source/Debris_Temperature_Profiles/Profile_Temperature_Records'/f"{p['profile_id']}.csv"
  assert sha(rawpath)==p['source_sha256']
  with rawpath.open(encoding='utf-8-sig') as f:raw=list(csv.DictReader(f))
  rawmap={r['Datetime'].replace(' ','T')+'Z':r for r in raw}
  probe={float(q['depth_m']):q for q in p['per_probe']};train0=dt.datetime.fromisoformat(p['split']['initial_state_utc'].rstrip('Z'))
  grouped={}
  for row in rows:
   z=float(row['depth_m']);t=dt.datetime.fromisoformat(row['utc'].rstrip('Z'));seconds=(t-train0).total_seconds()
   assert seconds%3600==0 and seconds>=86400
   expected='train' if seconds<14*86400 else 'test';assert row['role']==expected and seconds<=28*86400
   key=(p['profile_id'],row['utc'],z);assert key not in allids;allids.add(key)
   source=rawmap[row['utc']];columns=[c for c in source if c.startswith('Depth:')]
   lookup={float(c.split(':')[1].split('cm')[0])/100:c for c in columns}
   assert float(source[lookup[z]])==float(row['observed_C'])
   assert float(source[lookup[min(lookup)]])==float(row['upper_boundary_C'])
   assert float(source[lookup[max(lookup)]])==float(row['lower_boundary_C'])
   assert all(math.isfinite(float(row[m+'_C'])) for m in MODELS)
   grouped.setdefault((z,expected),[]).append(row)
  recomputed={}
  for z,qp in probe.items():
   for role,storedkey in [('train','training_metrics'),('test','test_metrics')]:
    rr=grouped[z,role];assert len(rr)==(312 if role=='train' else 337)
    for m in MODELS:
     mm=metric([float(r['observed_C']) for r in rr],[float(r[m+'_C']) for r in rr])
     for k,v in mm.items():diffs.append(abs(v-qp[storedkey][m][k]))
     if role=='test':recomputed.setdefault(m,[]).append(mm)
  profile={m:{k:math.fsum(q[k] for q in qq)/len(qq) for k in ['MAE_C','MSE_C','bias_C']} for m,qq in recomputed.items()}
  for m in MODELS:
   profile[m]['RMSE_C']=math.sqrt(profile[m]['MSE_C'])
   for k,v in profile[m].items():diffs.append(abs(v-p['test_metrics'][m][k]))
  assert p['alpha_best_grid_index']==min(range(len(p['alpha_grid_training_MSE_C2'])),key=lambda j:p['alpha_grid_training_MSE_C2'][j])
  train_mse=math.fsum(q['training_metrics']['diffusion']['MSE_C'] for q in p['per_probe'])/len(p['per_probe'])
  assert train_mse<=min(p['alpha_grid_training_MSE_C2'])+1e-10
  lam=min(p['ridge_inner_validation_MSE_C2'],key=p['ridge_inner_validation_MSE_C2'].get)
  assert float(lam)==p['ridge_penalty']
  profiles.append(dict(profile_id=p['profile_id'],metrics=profile))
 aggregate={}
 for m in MODELS:
  aggregate[m]={key:math.fsum(p['metrics'][m][k] for p in profiles)/len(profiles) for key,k in [('equal_profile_mean_MAE_C','MAE_C'),('equal_profile_mean_bias_C','bias_C')]}
  aggregate[m]['equal_profile_RMSE_C']=math.sqrt(math.fsum(p['metrics'][m]['MSE_C'] for p in profiles)/len(profiles))
  for k,v in aggregate[m].items():diffs.append(abs(v-result['aggregate_over_completed_profiles'][m][k]))
 assert len(profiles)==11 and len(allids)==17523 and max(diffs)<1e-10
 audit=dict(status='PASS',scope='saved-event arithmetic and source-truth membership; not independent physical validation',profiles=11,interior_series=27,unique_prediction_rows=len(allids),max_metric_difference=max(diffs),source_result_sha256=sha(ROOT/'results/FIELD_TEMPERATURE_RESULTS.json'),recomputed_aggregate=aggregate,checks=['raw source truth and endpoint matching','common hourly train/test timestamps','no duplicate prediction identities','all four models finite on identical rows','grid/penalty selection agrees with saved training objectives','equal-probe/profile weighting independently recomputed'])
 (ROOT/'results/ARITHMETIC_AND_MEMBERSHIP_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
if __name__=='__main__':run()
