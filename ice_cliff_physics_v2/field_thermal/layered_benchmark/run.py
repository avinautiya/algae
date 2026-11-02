import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import csv,json,hashlib,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent;FIELD=ROOT.parent
sys.path.insert(0,str(FIELD/'temperature_benchmark'));import run_benchmark as old
sys.path.insert(0,str(FIELD/'layered_solver'));from layered import PreparedLayered
sys.path.insert(0,str(FIELD/'solver'));from dirichlet_slab_temperature_v1_executed import PreparedSlab
IDS=['KDT_2021a','KDT_2021b','KDT_2021d','KDT_2022a','KDT_2022b','KDT_2022c','KDT_2023a','KDT_2023b','KDT_2023c']
MODELS=['layered','homogeneous','quasisteady','lag_filter','ridge','homogeneous_fv'];GRID=np.logspace(-8,-5,9)
def load(s):
 p=old.source_path(s['profile_id']);assert old.sha(p)==s['source_sha256']
 with p.open(encoding='utf-8-sig',newline='') as f:r=csv.DictReader(f);rows=list(r);cols=[x for x in r.fieldnames if x.startswith('Depth:')]
 n=42*288+1;rows=rows[:n];assert len(rows)==n
 times=np.array([(old.stamp(r['Datetime'])-old.stamp(rows[0]['Datetime'])).total_seconds() for r in rows]);assert np.all(np.diff(times)==300)
 z=np.array(s['all_depths_cm'])/100;v=np.array([[float(r[c]) for c in cols] for r in rows]);assert np.isfinite(v).all()
 return times,z,v,[r['Datetime'] for r in rows]
def main():
 out=ROOT/'results';assert not (out/'RESULTS.json').exists(),'choose new version; do not overwrite'
 freeze=json.loads((ROOT/'FREEZE.json').read_text())
 for name,h in freeze['files'].items():assert old.sha(FIELD/name)==h,name
 splits={s['profile_id']:s for s in json.loads((FIELD/'protocol/TEMPERATURE_SPLITS.json').read_text())};ridgeall=json.loads((FIELD/'temperature_benchmark/results/RIDGE_PARAMETERS.json').read_text())['models']
 results=[]
 for pid in IDS:
  start=time.perf_counter();s=splits[pid];ts,z,v,dates=load(s);train=(ts>=86400)&(ts<14*86400)&(np.arange(len(ts))%12==0);test=(ts>=28*86400)&(np.arange(len(ts))%12==0)
  ntrain=14*288+1;ps=PreparedLayered(ts[:ntrain],v[:ntrain,0],v[:ntrain,-1],z[1:-1],z,v[0],cells=40,max_gap_s=300)
  costs=[]
  for a in GRID:
   for b in GRID:
    p=ps.predict(float(a),float(b));costs.append(dict(alpha_top=float(a),alpha_bottom=float(b),training_MSE_C2=old.losses(v[:ntrain,1:-1][train[:ntrain]],p[train[:ntrain]])))
  best=min(costs,key=lambda x:x['training_MSE_C2']);a=best['alpha_top'];b=best['alpha_bottom'];params=json.loads((FIELD/'temperature_benchmark/results/profiles'/f'{pid}.json').read_text());alpha=params['alpha_m2_s']
  psfull=PreparedLayered(ts,v[:,0],v[:,-1],z[1:-1],z,v[0],cells=40,max_gap_s=300)
  pde=PreparedSlab(ts,v[:,0],v[:,-1],z[1:-1],z,v[0],modes=64,max_gap_s=300)
  pred={'layered':psfull.predict(a,b),'homogeneous':pde.predict(alpha),'quasisteady':pde.linear,'homogeneous_fv':psfull.predict(alpha,alpha)}
  pred['lag_filter']=np.column_stack([old.lagged(pde.linear[:,j:j+1],v[0,j+1:j+2],fp['tau_s'])[:,0] for j,fp in enumerate(params['filter_parameters'])])
  idx=np.flatnonzero(test);fit={k:np.array(x) for k,x in ridgeall[pid].items() if k in ['mean','scale','keep','coef','ymean']};fit['keep']=fit['keep'].astype(bool);rp=np.full_like(v[:,1:-1],np.nan);rp[idx]=old.ridge_predict(fit,old.ridge_features(v,idx));pred['ridge']=rp
  fine=PreparedLayered(ts,v[:,0],v[:,-1],z[1:-1],z,v[0],cells=80,max_gap_s=300)
  diffs={'layered':float(np.max(np.abs(fine.predict(a,b)[test]-pred['layered'][test]))),'homogeneous_fv':float(np.max(np.abs(fine.predict(alpha,alpha)[test]-pred['homogeneous_fv'][test])))}
  metrics={m:dict(MAE_C=float(np.mean([old.metric(v[test,j+1],pred[m][test,j])['MAE_C'] for j in range(len(z)-2)])),MSE_C2=float(np.mean((pred[m][test]-v[test,1:-1])**2))) for m in MODELS}
  for m in metrics:metrics[m]['RMSE_C']=metrics[m]['MSE_C2']**.5
  out.mkdir(exist_ok=True);dest=out/f'{pid}.csv'
  with dest.open('w',newline='') as f:
   w=csv.writer(f);w.writerow(['profile','utc','depth_m','observed_C']+MODELS)
   for i in idx:
    for j,depth in enumerate(z[1:-1]):w.writerow([pid,dates[i],depth,v[i,j+1]]+[pred[m][i,j] for m in MODELS])
  r=dict(profile_id=pid,training_objectives=costs,chosen=best,alpha_at_boundary=a in (GRID[0],GRID[-1]) or b in (GRID[0],GRID[-1]),prior_homogeneous_alpha=alpha,test_hourly_times=int(test.sum()),interior_probes=len(z)-2,metrics=metrics,numerical_max40_vs80_C=diffs,numerical_pass=max(diffs.values())<.02,runtime_s=time.perf_counter()-start,prediction_sha256=old.sha(dest))
  old.put(out/f'{pid}.json',r);results.append(r);print(pid,'complete',flush=True)
 agg={m:dict(MAE_C=float(np.mean([r['metrics'][m]['MAE_C'] for r in results])),RMSE_C=float(np.mean([r['metrics'][m]['MSE_C2'] for r in results])**.5)) for m in MODELS};wins={m:sum(r['metrics']['layered']['MAE_C']<r['metrics'][m]['MAE_C'] for r in results) for m in MODELS if m!='layered'}
 old.put(out/'RESULTS.json',dict(status='COMPLETE',profiles=results,aggregate=agg,wins=wins,numerical_gate=all(r['numerical_pass'] for r in results),declared_gate_pass=all(wins[m]>=8 and agg['layered']['MAE_C']<agg[m]['MAE_C'] for m in ['homogeneous','quasisteady','lag_filter','ridge']) and all(r['numerical_pass'] for r in results),scope='reserved temporal observed-boundary interior temperature test; no melt inference',freeze_sha256=old.sha(ROOT/'FREEZE.json')));print(json.dumps(agg,indent=2),flush=True)
if __name__=='__main__':main()
