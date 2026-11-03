"""Causal rolling conductance hypothesis, development only."""
import csv,datetime as dt,json,math,sys,os
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from package_paths import round3_root
R3=round3_root()
sys.path[:0]=[str(R3/'theory'),str(R3/'moisture'),str(R3)]
from weak_graph import fit
from anomaly_index import load_site,dose
import importlib.util
_spec=importlib.util.spec_from_file_location('prior_paths',R3/'paths.py');_prior=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_prior)
OLD=_prior.require_piramide();C=1855487.;L=334000.

def phase_deliver(heat_J_m2):
    """Conditional cold debt; positive release only after repaying negative heat."""
    q=np.asarray(heat_J_m2,float)
    if q.ndim!=1 or not np.isfinite(q).all():raise ValueError('finite signed heat series required')
    cold=0.;out=[]
    for energy in q:
        state=cold+energy
        released=max(state,0.);cold=min(state,0.)
        out.append(released)
    return np.array(out),cold

def run_site(site,days):
    f=load_site(site);labels=f['labels'];n=len(labels)
    rows={r['SOURCE_LOCAL_LABEL']:r for r in csv.DictReader((OLD/'data/decoded'/('debtemp_'+site+'.csv')).open())}
    T=np.array([[float(rows[t.isoformat()]['DT'+str(j)])for j in range(1,6)]+[0.]for t in labels])
    if not np.isfinite(T).all():raise ValueError('missing temperature')
    # Absolute source depths: covered DT1=.02m, remaining nominal depths;
    # stipulated ice plane20mm below deepest probe. Coordinates relative toDT1.
    z=np.array([0.,.105,.23,.355,.48,.50])if site=='Tower1'else np.array([0.,.09,.20,.31,.42,.44])
    t=np.arange(n)*1800.;train_stop=f['slices'][0].stop
    first=fit(z,t[:train_stop+1],T[:train_stop+1],capacity_J_m3K=C,window_s=21600.)
    G=np.array(first['conductances_W_m2K']);history=[dict(index=train_stop,conductivities=first['conductivities_W_mK'],rmse=first['balance_RMSE_W_m2'])]
    q=np.empty(n-1);q[:train_stop]=G[-1]*(T[:train_stop,-2]+T[1:train_stop+1,-2])/2*1800
    for a in range(train_stop,n-1,48):
        start=max(0,a-int(days*48));solution=fit(z,t[start:a+1],T[start:a+1],capacity_J_m3K=C,window_s=21600.)
        G=np.array(solution['conductances_W_m2K']);b=min(a+48,n-1)
        q[a:b]=G[-1]*(T[a:b,-2]+T[a+1:b+1,-2])/2*1800
        history.append(dict(index=a,conductivities=solution['conductivities_W_mK'],rmse=solution['balance_RMSE_W_m2'],bounds_active=solution['physical_bound_active']))
    release,cold=phase_deliver(q);raw=np.array([sum(release[s])/L for s in f['slices']])
    return f,raw,dict(capacity_nominal=C,signed_heat_J_m2=float(sum(q)),latent_release_J_m2=float(sum(release)),final_cold_debt_J_m2=cold,
        phase_ledger_residual_J_m2=float(sum(q)-sum(release)-cold),source_depths_m=z.tolist(),conductance_updates=history)

def run():
    d=json.loads((R3/'moisture/DEVELOPMENT_FIT.json').read_text());p=json.loads((R3/'moisture/frozen_tower2_final/PREDICTIONS.json').read_text());r=json.loads((R3/'moisture/frozen_tower2_final/RESULTS.json').read_text())
    alltrials=[]
    for site in ['Tower1','Tower2']:
        y0=d['train_mmwe']if site=='Tower1'else p['local_train_target']['observed_ice_ablation_mmwe']
        truth=np.array(d['development_mmwe']if site=='Tower1'else [x['observed_mmwe']for x in r['model_metrics']['air_pdd']['period_errors']])
        for days in [7,3,14]:
            f,raw,diagnosis=run_site(site,days)
            if raw[0]<=0:raise ValueError('no positive training heat, refuse melt calibration')
            gain=y0/raw[0];pred=raw[1:]*gain;air=y0*dose(f,'air_pdd')[1:]/dose(f,'air_pdd')[0];cal=y0*dose(f,'calendar_exp_air',-.0031)[1:]/dose(f,'calendar_exp_air',-.0031)[0]
            trial=dict(site=site,window_days=days,primary=days==7,predictions_mmwe=pred.tolist(),observed_mmwe=truth.tolist(),MAE_mmwe=float(np.mean(abs(pred-truth))),air_MAE_mmwe=float(np.mean(abs(air-truth))),calendar_MAE_mmwe=float(np.mean(abs(cal-truth))),
                heat_scale_from_TRAIN=gain,implied_scaled_capacity_J_m3K=C*gain,implied_capacity_within_declared_1_to_3MJ=bool(1e6<=C*gain<=3e6),raw_period_releases_mmwe=raw.tolist(),diagnosis=diagnosis)
            alltrials.append(trial);print(site,days,trial['MAE_mmwe'],trial['implied_scaled_capacity_J_m3K'],flush=True)
    with(HERE/'RESULTS.json').open('x')as out:json.dump(dict(scope='Already inspected DEVELOPMENT only',historical_novelty=False,field_superiority=False,trials=alltrials),out,indent=2,allow_nan=False);out.write('\n')
if __name__=='__main__':run()
