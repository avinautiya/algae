"""Independent saved-result/raw-data audit, fixed-parameter replay, no fitting."""
import datetime as dt,hashlib,importlib.util,json,math,statistics,sys
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent/'extension';OUT=ROOT/'results/frozen_point_diagnostic'
DATA=Path('/private/tmp/glacier-next-20261103/benchmark/CHANGRI_ALIGNED_AUDIT.json')
START=dt.datetime.fromisoformat('2015-12-01T00:00:00+05:45').timestamp()
CUT=dt.datetime.fromisoformat('2015-12-29T00:00:00+05:45').timestamp()
END=dt.datetime.fromisoformat('2016-01-12T00:00:00+05:45').timestamp()
REMOTE={'FREEZE.json':'446cbd08017f70560559757371e2f1e0206bf1df','plate_experiment.py':'7cb0f7a08ac4530e492dcab2fd85b777559fb426','dirichlet_slab.py':'a7b48233513a7025a368504b303168cfd3fd8ac3','PROTOCOL.md':'91c3c5a198f2cc317fa7014c4efa5f1aab737139'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def blob(p):
 b=p.read_bytes();return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def epoch(s):return dt.datetime.fromisoformat(s).timestamp()
def metrics(y,z):
 e=[b-a for a,b in zip(y,z)];n=len(e)
 if not n or len(y)!=len(z) or not all(math.isfinite(v) for v in e):raise AssertionError('invalid comparison')
 return dict(n=n,mae_W_m2=math.fsum(map(abs,e))/n,rmse_W_m2=math.sqrt(math.fsum(v*v for v in e)/n),bias_W_m2=math.fsum(e)/n)
def main():
 freeze=json.loads((ROOT/'FREEZE.json').read_text());result=json.loads((OUT/'RESULTS.json').read_text());pred=json.loads((OUT/'PREDICTIONS.json').read_text())
 assert result['preexecution_commit']=='2a0828dd2a4c8c42562a277eec8dfb2b991cf035'
 assert all(blob(ROOT/n)==v for n,v in REMOTE.items())
 assert all(sha(ROOT/n)==v for n,v in freeze['source_sha256'].items())
 assert sha(DATA)==freeze['data_sha256']
 raw=json.loads(DATA.read_text())['hourly'];raw=[r for r in raw if epoch(r['timestamp_utc'])<END]
 times=np.array([epoch(r['timestamp_utc']) for r in raw]);assert np.allclose(np.diff(times),3600,rtol=0,atol=1e-5)
 lookup={r['timestamp_utc']:r for r in raw};assert len(lookup)==len(raw)
 seen=set();rows_by_site={s:[] for s in [1,2]}
 for r in pred:
  key=(r['site'],r['timestamp_utc']);assert key not in seen;seen.add(key);assert CUT<=epoch(r['timestamp_utc'])<END
  actual=lookup[r['timestamp_utc']]['interface_flux_W_m2']['site'+str(r['site'])]
  assert actual==r['observed_raw_plate_W_m2'];rows_by_site[r['site']].append(r)
 # Replay only the frozen solver and saved parameters. No parameter is fitted.
 spec=importlib.util.spec_from_file_location('audited_frozen_slab',ROOT/'dirichlet_slab.py');module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
 audits=[];max_metric_diff=0.;max_replay_diff=0.
 for saved in result['results']:
  site=saved['site'];rr=rows_by_site[site];expected=[r['timestamp_utc'] for r in raw if CUT<=epoch(r['timestamp_utc'])<END and r['interface_flux_W_m2'] is not None and r['interface_flux_W_m2'].get('site'+str(site)) is not None]
  assert [r['timestamp_utc'] for r in rr]==expected and len(rr)==336
  allmetrics={}
  for name in rr[0]['predictions']:
   m=metrics([r['observed_raw_plate_W_m2'] for r in rr],[r['predictions'][name] for r in rr]);allmetrics[name]=m
   for key,value in m.items():max_metric_diff=max(max_metric_diff,abs(value-saved['scores'][name][key]))
  names=['T3 -1cm','T2 -5cm','T1 -10cm'] if site==1 else ['T5 -1cm','T4 -8cm'];depths=[.01,.05,.10] if site==1 else [.01,.08]
  temp=np.array([[r['probe_temperature_C'][n] for n in names] for r in raw]);targets=[.05] if site==1 else [.045]
  slab=module.PreparedSlab(times-times[0],temp[:,0],temp[:,-1],targets,depths,temp[0],modes=64,max_gap_s=3600.)
  slab128=module.PreparedSlab(times-times[0],temp[:,0],temp[:,-1],targets,depths,temp[0],modes=128,max_gap_s=3600.)
  coef=saved['fits']['diffusion']['raw_signal_coefficient_W_m_K'];alpha=saved['alpha_m2_s'];prediction=coef*slab.bottom_downward_gradient(alpha).gradient_K_m;fine=coef*slab128.bottom_downward_gradient(alpha).gradient_K_m
  features=np.empty((len(times),12));features[:]=np.nan
  for j,lag in enumerate([0,1,3,6,12,24]):
   features[lag:,2*j:2*j+2]=temp[:len(times)-lag if lag else len(times),:][:,[0,-1]]
  fit=saved['fits']['ridge'];ridge=np.column_stack([np.ones(len(times)),(features-np.array(fit['mean']))/np.array(fit['scale'])])@np.array(fit['coefficients'])
  train=(times>=START)&(times<CUT);training_q=[r['interface_flux_W_m2']['site'+str(site)] for i,r in enumerate(raw) if train[i] and r['interface_flux_W_m2'] is not None and r['interface_flux_W_m2'].get('site'+str(site)) is not None];median=statistics.median(training_q)
  ii={r['timestamp_utc']:i for i,r in enumerate(raw)};indices=[ii[r['timestamp_utc']] for r in rr]
  for r,i in zip(rr,indices):
   instant=saved['fits']['instantaneous']['raw_signal_coefficient_W_m_K']*(temp[i,0]-temp[i,-1])/(depths[-1]-depths[0])
   replay=dict(diffusion=prediction[i],instantaneous_gradient=instant,training_median=median,causal_boundary_ridge=ridge[i],diffusion_opposite_orientation_sensitivity=-prediction[i])
   for name,z in replay.items():max_replay_diff=max(max_replay_diff,abs(float(z)-r['predictions'][name]))
  numeric=float(np.max(abs(prediction[indices]-fine[indices])));assert abs(numeric-saved['numerics']['max64_vs128_signal_difference_W_m2'])<1e-10
  trainrms=math.sqrt(math.fsum(v*v for v in training_q)/len(training_q));limit=max(.5,.01*trainrms);assert abs(limit-saved['numerics']['limit_W_m2'])<1e-12;assert (numeric<=limit)==saved['numerics']['pass_']
  assert saved['alpha_m2_s']==min(result['alpha_training_search'])[1]
  audits.append(dict(site=site,rows=len(rr),scores=allmetrics,numerical_difference_W_m2=numeric,numerical_limit_W_m2=limit,numerical_gate_pass=numeric<=limit))
 mean=lambda name:math.fsum(a['scores'][name]['mae_W_m2'] for a in audits)/2
 pde=mean('diffusion');ridge=mean('causal_boundary_ridge');both=all(a['scores']['diffusion']['mae_W_m2']<a['scores']['causal_boundary_ridge']['mae_W_m2'] for a in audits)
 gate=both and pde<=.95*ridge and all(a['numerical_gate_pass'] for a in audits) and not any(r['fitted_alpha_is_grid_endpoint'] for r in result['results'])
 assert gate==result['frozen_superiority_rule_pass'];assert abs(pde-result['primary_equal_site_mae_W_m2'])<1e-12;assert abs(ridge-result['ridge_equal_site_mae_W_m2'])<1e-12
 assert max_metric_diff<1e-10 and max_replay_diff<1e-10 and len(seen)==672
 audit=dict(preexecution_commit=result['preexecution_commit'],remote_freeze_and_sources_git_blob_hashes_verified=True,data_sha256_verified=True,raw_target_rows_checked=len(seen),fixed_support_pass=True,fit_parameters_refitted=False,max_metric_recomputation_difference=max_metric_diff,max_fixed_parameter_prediction_replay_difference=max_replay_diff,primary_equal_site_diffusion_MAE_W_m2=pde,primary_equal_site_ridge_MAE_W_m2=ridge,diffusion_relative_MAE_increase_vs_ridge=pde/ridge-1,frozen_superiority_rule_pass=gate,sites=audits,scope='conditional representative-point raw plate signal only; not melt, native interval heat flux, novel physics or environmental superiority')
 (HERE/'POSTSCORE_AUDIT.json').write_text(json.dumps(audit,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in audit.items() if k!='sites'},indent=2))
if __name__=='__main__':main()
