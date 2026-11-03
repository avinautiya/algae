"""Conditional measured-boundary plate-signal experiment; not melt forecasting."""
import argparse,datetime as dt,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from dirichlet_slab import PreparedSlab
ROOT=Path(__file__).resolve().parent
START='2015-12-01T00:00:00+05:45';CUT='2015-12-29T00:00:00+05:45';END='2016-01-12T00:00:00+05:45'
ALPHAS=np.logspace(-8,-5,61)
def epoch(s):return dt.datetime.fromisoformat(s).timestamp()
def coefficient(x,y):
 x=np.asarray(x,float);y=np.asarray(y,float)
 if x.shape!=y.shape or not len(x) or not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('finite matching training vectors required')
 s=float(x@x)
 if s<1e-12:raise ValueError('unidentified flux amplitude: no gradient excitation')
 raw=float(x@y/s)
 if raw==0:raise ValueError('zero flux/gradient amplitude; orientation unidentified')
 return dict(raw_signal_coefficient_W_m_K=raw,effective_conductivity_abs_W_m_K=abs(raw),training_inferred_orientation=1 if raw>0 else -1,gradient_squared_sum=s)
def validate_hourly(times):
 a=np.asarray(times,float)
 if len(a)<2 or not np.isfinite(a).all() or not np.allclose(np.diff(a),3600.,rtol=0.,atol=1e-5):raise ValueError('exact hourly support required for causal hour lags')

def prepare(data,site,modes):
 names=['T3 -1cm','T2 -5cm','T1 -10cm'] if site==1 else ['T5 -1cm','T4 -8cm'];depths=[.01,.05,.10] if site==1 else [.01,.08]
 rows=[r for r in data['hourly'] if epoch(r['timestamp_utc'])<epoch(END)]
 # Fixed first source row initialization; never reset at training/test boundary.
 times=np.array([epoch(r['timestamp_utc']) for r in rows]);temp=np.array([[np.nan if r['probe_temperature_C'] is None or r['probe_temperature_C'].get(n) is None else r['probe_temperature_C'][n] for n in names] for r in rows])
 validate_hourly(times)
 if not np.isfinite(temp).all():raise ValueError('unresolved measured-boundary support')
 q=np.array([np.nan if r['interface_flux_W_m2'] is None or r['interface_flux_W_m2'].get('site'+str(site)) is None else r['interface_flux_W_m2']['site'+str(site)] for r in rows])
 target=[.05] if site==1 else [.045]
 slab=PreparedSlab(times-times[0],temp[:,0],temp[:,-1],target,depths,temp[0],modes=modes,max_gap_s=3600.)
 train=(times>=epoch(START))&(times<epoch(CUT));test=(times>=epoch(CUT))&(times<epoch(END))
 return dict(rows=rows,times=times,temp=temp,q=q,slab=slab,train=train,test=test,gradient=(temp[:,0]-temp[:,-1])/(depths[-1]-depths[0]))
def lagged(temp):
 lags=[0,1,3,6,12,24];features=[]
 for lag in lags:
  a=np.empty_like(temp[:,[0,-1]]);a[:lag]=np.nan
  if lag:a[lag:]=temp[:-lag][:,[0,-1]]
  else:a[:]=temp[:,[0,-1]]
  features.append(a)
 return np.column_stack(features)
def ridge_fit(X,y,mask,lam):
 xx=X[mask];yy=y[mask];mu=xx.mean(0);scale=xx.std(0);scale[scale<1e-9]=1
 z=np.column_stack([np.ones(len(xx)),(xx-mu)/scale]);pen=np.eye(z.shape[1])*lam;pen[0,0]=0
 b=np.linalg.lstsq(z.T@z+pen,z.T@yy,rcond=None)[0]
 pred=np.column_stack([np.ones(len(X)),(X-mu)/scale])@b
 return pred,dict(lambda_=lam,mean=mu.tolist(),scale=scale.tolist(),coefficients=b.tolist())
def models_for_site(p,alpha):
 q=p['q'];train=p['train']&np.isfinite(q);g=p['slab'].bottom_downward_gradient(alpha).gradient_K_m
 a=coefficient(g[train],q[train]);b=coefficient(p['gradient'][train],q[train]);pred=dict(diffusion=a['raw_signal_coefficient_W_m_K']*g,instantaneous_gradient=b['raw_signal_coefficient_W_m_K']*p['gradient'],training_median=np.full(len(q),np.median(q[train])))
 X=lagged(p['temp'][:,[0,-1]]);valid=np.isfinite(X).all(1);inner=train&(p['times']<epoch(CUT)-7*86400)&valid;validation=train&~inner&valid
 choices=[]
 for lam in [0.,.01,1.,100.]:
  z,_=ridge_fit(X,q,inner,lam);choices.append((float(np.mean(abs(z[validation]-q[validation]))),lam))
 selected=min(choices)[1];z,fit=ridge_fit(X,q,train&valid,selected);pred['causal_boundary_ridge']=z;pred['diffusion_opposite_orientation_sensitivity']=-pred['diffusion']
 return pred,dict(diffusion=a,instantaneous=b,ridge=fit,ridge_inner_validation=choices)
def scores(obs,pred,mask):
 y=obs[mask];z=pred[mask]
 if not len(y) or not np.isfinite(y).all() or not np.isfinite(z).all():raise ValueError('empty/nonfinite fixed test support')
 e=z-y;return dict(n=len(e),mae_W_m2=float(np.mean(abs(e))),rmse_W_m2=float(np.sqrt(np.mean(e*e))),bias_W_m2=float(np.mean(e)))
def run(data_path,out,commit,interpretation):
 if interpretation!='representative_point_diagnostic':raise ValueError('native interval-flux validation is blocked; explicit representative-point diagnostic required')
 if out.exists():raise ValueError('new output directory required')
 freeze=json.loads((ROOT/'FREEZE.json').read_text())
 if any(hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h for n,h in freeze['source_sha256'].items()):raise ValueError('source changed after freeze')
 if hashlib.sha256(data_path.read_bytes()).hexdigest()!=freeze['data_sha256']:raise ValueError('source data changed')
 if len(commit)!=40 or any(c not in '0123456789abcdef' for c in commit):raise ValueError('pre-execution published commit required')
 data=json.loads(data_path.read_text());p1=prepare(data,1,64);tt=p1['temp'][:,1];train=p1['train'];search=[]
 for alpha in ALPHAS:search.append((float(np.mean(abs(p1['slab'].predict(float(alpha))[:,0][train]-tt[train]))),float(alpha)))
 alpha=min(search)[1];results=[];table=[]
 for site in [1,2]:
  p=p1 if site==1 else prepare(data,site,64);pred,fit=models_for_site(p,alpha);p128=prepare(data,site,128);g128=p128['slab'].bottom_downward_gradient(alpha).gradient_K_m;fine=fit['diffusion']['raw_signal_coefficient_W_m_K']*g128
  mask=p['test']&np.isfinite(p['q']);numeric_error=float(np.max(abs(pred['diffusion'][mask]-fine[mask])));limit=max(.5,.01*float(np.sqrt(np.mean(p['q'][p['train']&np.isfinite(p['q'])]**2))))
  results.append(dict(site=site,alpha_m2_s=alpha,alpha_source='site1 training-only interior-temperature fit',fits=fit,scores={name:scores(p['q'],value,mask) for name,value in pred.items()},numerics=dict(max64_vs128_signal_difference_W_m2=numeric_error,limit_W_m2=limit,pass_=numeric_error<=limit),fitted_alpha_is_grid_endpoint=alpha in [float(ALPHAS[0]),float(ALPHAS[-1])]))
  for i in np.flatnonzero(mask):table.append(dict(site=site,timestamp_utc=p['rows'][i]['timestamp_utc'],observed_raw_plate_W_m2=float(p['q'][i]),predictions={k:float(v[i]) for k,v in pred.items()}))
 primary_mae=float(np.mean([r['scores']['diffusion']['mae_W_m2'] for r in results]));ridge_mae=float(np.mean([r['scores']['causal_boundary_ridge']['mae_W_m2'] for r in results]));primary_pass=all(r['scores']['diffusion']['mae_W_m2']<r['scores']['causal_boundary_ridge']['mae_W_m2'] and r['numerics']['pass_'] and not r['fitted_alpha_is_grid_endpoint'] for r in results) and ridge_mae>0 and primary_mae<=.95*ridge_mae
 out.mkdir(parents=True);(out/'RESULTS.json').write_text(json.dumps(dict(preexecution_commit=commit,results=results,alpha_training_search=search,primary_equal_site_mae_W_m2=primary_mae,ridge_equal_site_mae_W_m2=ridge_mae,frozen_superiority_rule_pass=primary_pass,interpretation='Conditional prediction of archived raw plate signal at representative timestamps. Native heat-flux averaging/sign/placement not independently resolved. NOT melt validation or a weather forecast.'),indent=2,allow_nan=False)+'\n');(out/'PREDICTIONS.json').write_text(json.dumps(table,indent=2,allow_nan=False)+'\n')
 print(json.dumps(dict(sites=len(results),predictions=len(table),numerical_gates=[r['numerics']['pass_'] for r in results])))
def freeze(data_path):
 p=ROOT/'FREEZE.json'
 if p.exists():raise ValueError('do not overwrite freeze')
 source={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['plate_experiment.py','dirichlet_slab.py','PROTOCOL.md']}
 d=json.loads(data_path.read_text());summary={}
 for site in [1,2]:
  p0=prepare(d,site,64);summary[str(site)]=dict(train_rows=int(np.sum(p0['train']&np.isfinite(p0['q']))),test_rows=int(np.sum(p0['test']&np.isfinite(p0['q']))))
 p.write_text(json.dumps(dict(source_sha256=source,data_sha256=hashlib.sha256(data_path.read_bytes()).hexdigest(),support_counts_only=summary,heldout_scores_not_computed=True,first14_source_days_training_diagnostic_already_seen=True),indent=2)+'\n');print('freeze written; no held-out score computed')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('data',type=Path);p.add_argument('--freeze',action='store_true');p.add_argument('--out',type=Path);p.add_argument('--commit');p.add_argument('--interpretation',choices=['representative_point_diagnostic']);a=p.parse_args()
 if a.freeze:freeze(a.data)
 elif a.out and a.commit:run(a.data,a.out,a.commit,a.interpretation)
 else:p.error('use --freeze or --out NEW --commit SHA')
