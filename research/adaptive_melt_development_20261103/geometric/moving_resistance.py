"""Classical quasi-steady moving-resistance development falsifier.

Not novel physics, not a confirmed account of debris settlement or field geometry.
No unknown target files are opened. Only released prior result JSONs are used.
"""
import json,math,sys,os
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from package_paths import round3_root
R3=round3_root()
sys.path.insert(0,str(R3/'moisture'))
from anomaly_index import load_site,dose
RHO=917.;LATENT=334000.

def incremental_melt(doses,train_mass_mmwe,depth_m,rho=RHO):
    f=np.asarray(doses,float)
    if f.ndim!=1 or len(f)<2 or not np.isfinite(f).all()or(f<0).any()or f[0]<=0:raise ValueError('finite nonnegative dose with positive TRAIN required')
    if not all(math.isfinite(x)and x>0 for x in [train_mass_mmwe,depth_m,rho]):raise ValueError('positive finite geometry and TRAIN mass required')
    # kg/m2 == mm water equivalent. Ice recession y/rho has units metres.
    recession=train_mass_mmwe/rho
    B=(2*depth_m*recession+recession**2)/f[0]
    cumulative=np.concatenate(([0.],np.cumsum(f)))
    H=np.sqrt(depth_m**2+B*cumulative)
    # Stable sqrt difference avoids cancellation when f is small.
    melt=rho*B*f/(H[1:]+H[:-1])
    k=RHO*LATENT*B/2
    return dict(period_melt_mmwe=melt.tolist(),depth_end_m=H.tolist(),inferred_k_W_mK=k,
                squared_depth_ledger_residual=float(np.max(np.abs(np.diff(H**2)-B*f))))

def run():
    dev=json.loads((R3/'moisture/DEVELOPMENT_FIT.json').read_text())
    pred=json.loads((R3/'moisture/frozen_tower2_final/PREDICTIONS.json').read_text())
    res=json.loads((R3/'moisture/frozen_tower2_final/RESULTS.json').read_text())
    result=dict(scope='Already inspected two-stake DEVELOPMENT only; classical moving thermal resistance hypothesis, not validated or historical novelty',sites={})
    for site,H in [('Tower1',.5),('Tower2',.45)]:
        f=load_site(site)
        if site=='Tower1':y0=dev['train_mmwe'];y=np.array(dev['development_mmwe'])
        else:y0=pred['local_train_target']['observed_ice_ablation_mmwe'];y=np.array([r['observed_mmwe']for r in res['model_metrics']['air_pdd']['period_errors']])
        trials=[]
        for boundary,method in [('primary_covered_surface','buried_pdd'),('secondary_empirical_air','air_pdd')]:
            q=dose(f,method);r=incremental_melt(q,y0,H)
            z=np.array(r['period_melt_mmwe'][1:]);fixed=y0*q[1:]/q[0]
            trials.append(dict(boundary=boundary,nominal_depth_m=H,predicted_mmwe=z.tolist(),observed_mmwe=y.tolist(),MAE_mmwe=float(np.mean(abs(z-y))),
                fixed_resistance_same_boundary_MAE_mmwe=float(np.mean(abs(fixed-y))),**r))
        a=y0*dose(f,'air_pdd')[1:]/dose(f,'air_pdd')[0]
        c=y0*dose(f,'calendar_exp_air',-.0031)[1:]/dose(f,'calendar_exp_air',-.0031)[0]
        result['sites'][site]=dict(trials=trials,air_MAE_mmwe=float(np.mean(abs(a-y))),calendar_MAE_mmwe=float(np.mean(abs(c-y))),epochs=f['epochs'])
    out=HERE/'RESULTS.json'
    with out.open('x')as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':run()
