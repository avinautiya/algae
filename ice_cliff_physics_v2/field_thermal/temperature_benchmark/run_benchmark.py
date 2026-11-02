"""Execute the frozen observed-boundary field-temperature test without test tuning."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,csv,datetime as dt,hashlib,json,math,platform,re,resource,sys,time,traceback
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import lfilter

ROOT=Path(__file__).resolve().parent
FIELD=ROOT.parent
PROTOCOL=FIELD/'protocol'
sys.path.insert(0,str(FIELD/'solver'))
from dirichlet_slab import PreparedSlab
MODELS=['diffusion','quasisteady','lag_filter','ridge']
ALPHAS=np.logspace(-8,-5,61)
TAUS=np.r_[0.,np.geomspace(60.,7*86400.,61)]
PENALTIES=[1e-6,1e-3,1.,1e3]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def put(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def stamp(s):return dt.datetime.fromisoformat(s.rstrip('Z')).replace(tzinfo=dt.timezone.utc)

def source_path(profile):
    return FIELD/'source/Debris_Temperature_Profiles/Profile_Temperature_Records'/f'{profile}.csv'

def load_profile(split):
    p=source_path(split['profile_id'])
    if sha(p)!=split['source_sha256']:raise ValueError('source hash changed '+split['profile_id'])
    with p.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f);headers=reader.fieldnames
        cols=[h for h in headers if re.fullmatch(r'Depth:\s*[\d.]+\s*cm',h)]
        depths=np.array([float(re.search(r'([\d.]+)\s*cm',h).group(1))/100 for h in cols])
        rows=list(reader)
    if not np.allclose(depths,np.array(split['all_depths_cm'])/100,rtol=0,atol=1e-12):raise ValueError('depth metadata changed')
    dates=[stamp(r['Datetime']) for r in rows]
    first=stamp(split['initial_state_utc']);end=stamp(split['test_end_utc'])
    if dates[0]!=first:raise ValueError('first input timestamp changed')
    use=[i for i,t in enumerate(dates) if t<=end]
    if not use or dates[use[-1]]!=end:raise ValueError('incomplete frozen horizon')
    dates=[dates[i] for i in use]
    times=np.array([(t-first).total_seconds() for t in dates])
    if len(times)!=8065 or not np.all(np.diff(times)==300):raise ValueError('native cadence/support differs from frozen audit')
    values=np.array([[float(rows[i][c]) for c in cols] for i in use])
    if not np.isfinite(values).all():raise ValueError('nonfinite sensor value; no silent filtering')
    for i in [use[0],use[-1]]:
        if (stamp(rows[i]['Datetime'])-stamp(rows[i]['Datetime_ADT'])).total_seconds()!=8*3600:raise ValueError('source clock relation changed')
    if not np.allclose(depths[1:-1],np.array(split['target_depths_cm'])/100):raise ValueError('interior membership changed')
    return times,dates,depths,values

def masks(times):
    hourly=np.arange(len(times))%12==0
    train=hourly&(times>=86400)&(times<14*86400)
    test=hourly&(times>=14*86400)&(times<=28*86400)
    return train,test

def losses(y,p):return float(np.mean(np.mean((p-y)**2,axis=0)))

def lagged(eq,initial,tau,delta=300.):
    if tau==0:return eq.copy()
    x=delta/tau;decay=math.exp(-x);phi=-math.expm1(-x)/x
    forcing=(phi-decay)*eq[:-1]+(1-phi)*eq[1:]
    out=np.empty_like(eq);out[0]=initial
    out[1:]=lfilter([1.],[1.,-decay],forcing,axis=0,zi=(decay*initial)[None,:])[0]
    return out

def ridge_fit(X,Y,penalty):
    mean=X.mean(axis=0);scale=X.std(axis=0);keep=scale>=1e-12
    if not keep.any():raise ValueError('no variable boundary features')
    Z=(X[:,keep]-mean[keep])/scale[keep];ym=Y.mean(axis=0)
    coef=np.linalg.solve(Z.T@Z+penalty*np.eye(Z.shape[1]),Z.T@(Y-ym))
    return dict(mean=mean,scale=scale,keep=keep,coef=coef,ymean=ym)
def ridge_predict(fit,X):
    return ((X[:,fit['keep']]-fit['mean'][fit['keep']])/fit['scale'][fit['keep']])@fit['coef']+fit['ymean']
def ridge_features(values,idx):
    return np.column_stack([values[idx-lag*12,b] for lag in [0,1,3,6,12,24] for b in [0,-1]])

def metric(y,p):
    e=p-y;a=np.abs(e)
    return dict(n=int(len(y)),MAE_C=float(a.mean()),RMSE_C=float(np.sqrt(np.mean(e*e))),MSE_C=float(np.mean(e*e)),bias_C=float(e.mean()),
                fraction_within_0_2C=float(np.mean(a<=.2)),fraction_within_0_5C=float(np.mean(a<=.5)),fraction_within_1C=float(np.mean(a<=1.)))

def execute_one(split,out):
    start=time.perf_counter();times,dates,depths,v=load_profile(split);truth=v[:,1:-1]
    train,test=masks(times);idx=np.flatnonzero(train|test)
    prepared=PreparedSlab(times,v[:,0],v[:,-1],depths[1:-1],depths,v[0],modes=64,max_gap_s=300.)
    objectives=[]
    for alpha in ALPHAS:
        pred=prepared.predict(float(alpha));objectives.append(losses(truth[train],pred[train]))
    best=int(np.argmin(objectives));alpha=float(ALPHAS[best]);refinement=None
    if 0<best<len(ALPHAS)-1:
        def fun(logalpha):return losses(truth[train],prepared.predict(math.exp(logalpha))[train])
        opt=minimize_scalar(fun,bounds=(math.log(ALPHAS[best-1]),math.log(ALPHAS[best+1])),method='bounded',options={'xatol':.005})
        if opt.success and float(opt.fun)<=objectives[best]:alpha=math.exp(float(opt.x))
        refinement=dict(success=bool(opt.success),training_MSE_C2=float(opt.fun),objective_calls=int(opt.nfev))
    predicted={'diffusion':prepared.predict(alpha),'quasisteady':prepared.linear.copy()}
    filt=np.empty_like(truth);filter_params=[]
    for j in range(truth.shape[1]):
        costs=[]
        for tau in TAUS:
            candidate=lagged(prepared.linear[:,j:j+1],truth[0,j:j+1],float(tau))[:,0]
            costs.append(float(np.mean((candidate[train]-truth[train,j])**2)))
        k=int(np.argmin(costs));tau=float(TAUS[k]);filter_params.append(dict(depth_m=float(depths[j+1]),tau_s=tau,training_MSE_C2=costs[k],at_grid_edge=bool(k==len(TAUS)-1)))
        filt[:,j]=lagged(prepared.linear[:,j:j+1],truth[0,j:j+1],tau)[:,0]
    predicted['lag_filter']=filt
    X=ridge_features(v,idx);Y=truth[idx];inner_train=times[idx]<10*86400;inner_val=(times[idx]>=10*86400)&(times[idx]<14*86400)
    cv=[]
    for lam in PENALTIES:
        fit=ridge_fit(X[inner_train],Y[inner_train],lam)
        cv.append(losses(Y[inner_val],ridge_predict(fit,X[inner_val])))
    lam=PENALTIES[int(np.argmin(cv))]
    final=ridge_fit(X[train[idx]],Y[train[idx]],lam)
    rr=np.full_like(truth,np.nan);rr[idx]=ridge_predict(final,X);predicted['ridge']=rr
    if not all(np.isfinite(predicted[m][idx]).all() for m in MODELS):raise ValueError('nonfinite prediction on common rows')
    numerical={}
    for n in [32,128]:
        p=PreparedSlab(times,v[:,0],v[:,-1],depths[1:-1],depths,v[0],modes=n,max_gap_s=300.).predict(alpha)
        numerical[f'max_{n}_vs64_C_on_common_rows']=float(np.max(np.abs(p[idx]-predicted['diffusion'][idx])))
    numerical['passes_0_02C_64_vs128']=numerical['max_128_vs64_C_on_common_rows']<.02
    perprobe=[]
    for j,z in enumerate(depths[1:-1]):
        perprobe.append(dict(depth_m=float(z),negative_truth_test_count=int(np.sum(truth[test,j]<0)),test_metrics={m:metric(truth[test,j],predicted[m][test,j]) for m in MODELS},training_metrics={m:metric(truth[train,j],predicted[m][train,j]) for m in MODELS}))
    profile={m:dict(MAE_C=float(np.mean([p['test_metrics'][m]['MAE_C'] for p in perprobe])),MSE_C=float(np.mean([p['test_metrics'][m]['MSE_C'] for p in perprobe])),bias_C=float(np.mean([p['test_metrics'][m]['bias_C'] for p in perprobe]))) for m in MODELS}
    for m in MODELS:profile[m]['RMSE_C']=math.sqrt(profile[m]['MSE_C'])
    dest=out/'predictions'/f"{split['profile_id']}.csv";dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.with_suffix('.tmp').open('w',newline='') as f:
        fields=['profile_id','utc','depth_m','role','observed_C','upper_boundary_C','lower_boundary_C']+[m+'_C' for m in MODELS]
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for i in idx:
            for j,z in enumerate(depths[1:-1]):
                row=dict(profile_id=split['profile_id'],utc=dates[i].isoformat().replace('+00:00','Z'),depth_m=float(z),role='train' if train[i] else 'test',observed_C=float(truth[i,j]),upper_boundary_C=float(v[i,0]),lower_boundary_C=float(v[i,-1]))
                row.update({m+'_C':float(predicted[m][i,j]) for m in MODELS});writer.writerow(row)
    dest.with_suffix('.tmp').replace(dest)
    daily=[]
    for day in range(14):
        daymask=test&(times>=14*86400+day*86400)&(times<14*86400+(day+1)*86400)
        daily.append(dict(test_elapsed_day=day,source_hourly_times=int(daymask.sum()),metrics={m:metric(truth[daymask].ravel(),predicted[m][daymask].ravel()) for m in MODELS}))
    return dict(profile_id=split['profile_id'],source_sha256=split['source_sha256'],source_sensor_depths_m=depths.tolist(),split=split,
        initialization='first source row only; continuous evolution; no split/test reset',train_hourly_times=int(train.sum()),test_hourly_times=int(test.sum()),interior_probes=len(perprobe),
        alpha_m2_s=alpha,alpha_best_grid_index=best,alpha_grid=ALPHAS.tolist(),alpha_grid_training_MSE_C2=objectives,alpha_refinement=refinement,alpha_at_grid_edge=bool(best in [0,len(ALPHAS)-1]),
        initial_modal_error_C=prepared.initial_modal_error,filter_parameters=filter_params,ridge_penalty=lam,ridge_inner_validation_MSE_C2=dict(zip(map(str,PENALTIES),cv)),
        ridge_fitted_feature_count=int(final['keep'].sum()),ridge_standardization='training only; intercept unpenalized',per_probe=perprobe,test_metrics=profile,daily_blocks=daily,
        numerical=numerical,prediction_csv=dest.relative_to(out).as_posix(),prediction_sha256=sha(dest),runtime_s=time.perf_counter()-start)

def prepare_freeze(out):
    frozen=json.loads((PROTOCOL/'TEMPERATURE_PROTOCOL_FREEZE.json').read_text())
    for name,key in [('TEMPERATURE_VALIDATION_SPEC.md','spec_sha256'),('TEMPERATURE_SPLITS.json','splits_sha256')]:
        if sha(PROTOCOL/name)!=frozen[key]:raise ValueError('protocol input changed '+name)
    allsplits=json.loads((PROTOCOL/'TEMPERATURE_SPLITS.json').read_text());splits=[s for s in allsplits if s['primary_eligible']]
    if len(splits)!=11 or sum(len(s['target_depths_cm']) for s in splits)!=27:raise ValueError('frozen cohort changed')
    audit=json.loads((FIELD/'PROFILE_FEASIBILITY.json').read_text());lookup={p['profile']:p for p in audit['profiles']}
    checks={}
    for s in splits:
        a=lookup[s['profile_id']]
        checks[s['profile_id']]=dict(source_sha_identical=a['source_sha256']==s['source_sha256'],depths_identical=a['depth_columns_cm']==s['all_depths_cm'],raw_bytes_match=sha(source_path(s['profile_id']))==s['source_sha256'])
    if not all(all(c.values()) for c in checks.values()):raise ValueError('scientific inputs changed after metadata audit update')
    put(out/'PREFIT_AUDIT_UPDATE.json',dict(original_freeze_unchanged=True,frozen_audit_sha256=frozen['feasibility_sha256'],current_audit_sha256=sha(FIELD/'PROFILE_FEASIBILITY.json'),scientific_input_checks=checks,interpretation='audit metadata was updated; frozen raw data, depths, splits and scientific model specification are unchanged; no field scores computed before this record'))
    put(out/'EXECUTION_FREEZE.json',dict(protocol_freeze=frozen,runner_sha256=sha(__file__),solver_sha256=sha(FIELD/'solver/dirichlet_slab.py'),source_checks=checks,
        models=MODELS,alpha_grid=ALPHAS.tolist(),tau_grid=TAUS.tolist(),ridge_penalties=PENALTIES,modes=64,precision_modes=[32,128],fit_labels='training objective only',environment=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),threads={k:os.environ[k] for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']})))
    return splits

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=ROOT/'results');a=p.parse_args();out=a.out
    if (out/'FIELD_TEMPERATURE_RESULTS.json').exists():raise ValueError('output already contains final scores; choose a versioned directory')
    splits=prepare_freeze(out)
    # Fixed-parameter timing pilot: no objective/accuracy or fitted parameter output.
    ts,dates,depths,v=load_profile(splits[0]);pilot=PreparedSlab(ts,v[:,0],v[:,-1],depths[1:-1],depths,v[0],modes=64,max_gap_s=300.)
    t=time.perf_counter();pilot.predict(1e-6);elapsed=time.perf_counter()-t
    put(out/'TIMING_PILOT.json',dict(runtime_s=elapsed,native_times=len(ts),target_depths=len(depths)-2,modes=64,accuracy_not_computed=True,parameter_not_fit=True,peak_memory_ru_maxrss_native_units=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    print('Timing pilot completed; runtime_s='+str(elapsed),flush=True)
    done=[];failed=[]
    for s in splits:
        try:
            r=execute_one(s,out);done.append(r);put(out/'profiles'/f"{s['profile_id']}.json",r)
            print('Completed '+s['profile_id']+' in '+str(round(r['runtime_s'],3))+'s',flush=True)
        except Exception as exc:
            failed.append(dict(profile_id=s['profile_id'],error=str(exc),traceback=traceback.format_exc()));put(out/'FAILURES.json',failed)
            print('FAILED '+s['profile_id']+': '+str(exc),flush=True)
    aggregate={}
    if done:
        for m in MODELS:
            aggregate[m]=dict(equal_profile_mean_MAE_C=float(np.mean([r['test_metrics'][m]['MAE_C'] for r in done])),equal_profile_RMSE_C=math.sqrt(float(np.mean([r['test_metrics'][m]['MSE_C'] for r in done]))),equal_profile_mean_bias_C=float(np.mean([r['test_metrics'][m]['bias_C'] for r in done])))
    complete=len(done)==11 and not failed;numerical=complete and all(r['numerical']['passes_0_02C_64_vs128'] for r in done)
    wins={m:sum(r['test_metrics']['diffusion']['MAE_C']<r['test_metrics'][m]['MAE_C'] for r in done) for m in MODELS if m!='diffusion'}
    result=dict(scope='observed-boundary retrospective subcolumn temperature validation; not melt or wet-film validation',status='COMPLETE' if complete else 'INCOMPLETE',numerical_precision_gate_passes=numerical,
        declared_profiles=11,completed_profiles=len(done),interior_series=sum(r['interior_probes'] for r in done),test_predictions=sum(r['interior_probes']*r['test_hourly_times'] for r in done),profile_results=done,failures=failed,
        aggregate_over_completed_profiles=aggregate,diffusion_profiles_better_than_baseline=wins,
        consistent_diffusion_superiority=bool(complete and numerical and all(x==11 for x in wins.values())),
        exclusions=[s for s in json.loads((PROTOCOL/'TEMPERATURE_SPLITS.json').read_text()) if not s['primary_eligible']],
        interpretation=['one glacier; profiles/days/probes are dependent','first24h initialization/history;13day training objective then14day test','negative temperatures retained, no freezing phase model inferred','k and Cv remain unidentifiable in scale; no flux or melt accuracy follows','later measured endpoints are imposed forcing, not issued forecasts','source depth and sensor representativeness assumptions retained'])
    put(out/'FIELD_TEMPERATURE_RESULTS.json',result)
    print(json.dumps(dict(status=result['status'],numerical_gate=numerical,aggregate=aggregate,wins=wins),indent=2),flush=True)
    if not complete:raise SystemExit(1)
if __name__=='__main__':main()
