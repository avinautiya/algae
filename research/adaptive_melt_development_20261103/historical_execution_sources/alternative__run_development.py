"""Finite4timescale DEVELOPMENT campaign on already released T1/T2 only."""
import argparse,json,math,hashlib,sys
from pathlib import Path
import numpy as np,pandas as pd
from radiation_reservoir import *
R3=Path('/private/tmp/glacier-round3-20261103');R2=Path('/private/tmp/glacier-round2-20261103')
sys.path.insert(0,str(R3/'moisture'))
from anomaly_index import load_site,dose,prediction_ratios

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def simulate(tower,tau_days,dt_seconds=1800):
    temp=pd.read_csv(R2/f'data/decoded/debtemp_{tower}.csv',index_col=0,parse_dates=True)
    weather=pd.read_csv(R3/'data/moisture/PIR22-23_60-kd_ERA5.csv',index_col=0,parse_dates=True)
    end=pd.Timestamp('2023-02-10T03:00');temp=temp.loc[:end]
    origin=temp.index[0];labels=pd.date_range(origin,end,freq=f'{dt_seconds}s')
    assert labels[-1]==end
    tn=temp.index.as_unit('ns').asi8.astype(float);wn=weather.index.as_unit('ns').asi8.astype(float)
    mid=labels[:-1]+pd.Timedelta(seconds=dt_seconds/2);mn=mid.as_unit('ns').asi8.astype(float)
    assert mn.min()>=wn.min() and mn.max()<=wn.max()
    assert np.isfinite(weather[['SWIN','LWIN']]).all().all()
    assert np.max(np.diff(wn))/1e9<=3600
    sw=np.interp(mn,wn,weather.SWIN);lw=np.interp(mn,wn,weather.LWIN);Ts=np.interp(mn,tn,temp.DT1)
    Q=np.array([radiation_proxy(a,b,c) for a,b,c in zip(sw,lw,Ts)])
    if tower=='Tower1':d=.5;C=1855487.1033697424;z=np.array([0,.125,.25,.375,.5])
    else:d=.45;C=1805830.3312199248;z=np.array([0,.11,.22,.33,.44,.45])
    initial_values=temp.iloc[0].to_numpy();initial_values=np.r_[initial_values,initial_values[-1]] if len(z)==6 else initial_values
    E0=C*float(np.trapezoid(initial_values,z));E=E0;tau=tau_days*86400*(d/.5)**2
    delivered=[];energies=[E];residual=[]
    for q in Q:
        new,g=step(E,float(q),dt_seconds,tau);residual.append(q*dt_seconds-(new-E)-g);delivered.append(g);E=new;energies.append(E)
    f=load_site(tower);epochs=pd.to_datetime(f['epochs']);out=[]
    for a,b in zip(epochs[:-1],epochs[1:]):
        mask=(labels[:-1]>=a)&(labels[1:]<=b);assert int(mask.sum())*dt_seconds==(b-a).total_seconds()
        out.append(math.fsum(np.asarray(delivered)[mask])/LATENT)
    trace=pd.DataFrame({'SOURCE_LOCAL_LABEL':labels[:-1].astype(str),'net_radiation_proxy_W_m2':Q,'SWIN_W_m2':sw,'LWIN_W_m2':lw,'buried_temperature_C_proxy':Ts,'effective_E_start_J_m2':energies[:-1],'delivered_proxy_J_m2':delivered,'ledger_residual_J_m2':residual})
    return np.asarray(out),trace,{'initial_E_J_m2':E0,'final_E_J_m2':E,'tau_site_seconds':tau,'initial_C_is_training_proxy_not_measured':C,
     'all_steps_max_abs_ledger_J_m2':float(np.max(np.abs(residual))),'whole_run_ledger_residual_J_m2':float(math.fsum(Q*dt_seconds)-(E-E0)-math.fsum(delivered)),
     'negative_net_radiation_proxy_steps':int((Q<0).sum()),'negative_effective_storage_steps':int((np.array(energies[:-1])<0).sum()),
     'interpretation':'Effectiveenergyreservoirnumericalclosure only; H/LE/snow/trueLWsurfaceunknown; effectiveEnegativeisnotmeasuredicecoldcontent'}

def run(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    dev=json.loads((R3/'moisture/DEVELOPMENT_FIT.json').read_text());t2pred=json.loads((R3/'moisture/frozen_tower2_final/PREDICTIONS.json').read_text());t2=json.loads((R3/'moisture/frozen_tower2_final/RESULTS.json').read_text())
    truth={'Tower1':{'train':dev['train_mmwe'],'development':np.array(dev['development_mmwe'])},'Tower2':{'train':float(t2pred['local_train_target']['observed_ice_ablation_mmwe']),'development':np.array([x['observed_mmwe'] for x in t2['model_metrics']['air_pdd']['period_errors']])}}
    trials=[];fixed_baselines={};calendar_pool=[]
    for tower in ['Tower1','Tower2']:
        f=load_site(tower);Y=truth[tower]['development'];Y0=truth[tower]['train']
        fixed_baselines[tower]={}
        for method,p in [('air_pdd',0.),('calendar_exp_air',dev['fits']['calendar_exp_air']['selected']['parameter'])]:
            pred=Y0*prediction_ratios(f,method,p);fixed_baselines[tower][method]={'parameter':p,'predictions_mmwe':pred.tolist(),'DEVELOPMENT_MAE_mmwe':float(np.mean(abs(pred-Y))),'signed_errors_mmwe':(pred-Y).tolist()}
        for gamma in np.linspace(-.04,.04,801):calendar_pool.append({'tower':tower,'parameter':float(gamma),'predictions':Y0*prediction_ratios(f,'calendar_exp_air',gamma)})
        for tau in TIMESCALE_DAYS:
            h,trace,diag=simulate(tower,tau);row={'tower':tower,'timescale_reference_days':tau,'potential_delivered_mmwe_all5':h.tolist(),'diagnostics':diag}
            if h[0]<=0:row.update({'supported':False,'reason':'NonpositiveTRAINpotential'});trials.append(row);continue
            gain=Y0/h[0];row['gain_fraction_TRAIN']=gain
            if not 0<=gain<=1:row.update({'supported':False,'reason':'TRAINpartitiongainoutside[0,1]'});trials.append(row);continue
            pred=gain*h[1:];row.update({'supported':True,'predictions_mmwe':pred.tolist(),'DEVELOPMENT_MAE_mmwe':float(np.mean(abs(pred-Y))),'signed_errors_mmwe':(pred-Y).tolist(),'calibrated_latent_mmwe_all5':(gain*h).tolist(),'unresolved_partition_fraction':1-gain})
            refined,_,_=simulate(tower,tau,900);rp=Y0*refined[1:]/refined[0];row['dt900_prediction_max_abs_difference_mmwe']=float(np.max(abs(rp-pred)))
            trace.to_csv(output/f'{tower}_tau{tau:g}_trace.csv',index=False);trials.append(row)
    supported=[x for x in trials if x['supported']];valid_tau=[tau for tau in TIMESCALE_DAYS if all(any(x['tower']==s and x['timescale_reference_days']==tau for x in supported) for s in ['Tower1','Tower2'])]
    if not valid_tau:raise ValueError('No supported familymember; preservedfailuretrialsneeded')
    cases=[(s,k) for s in ['Tower1','Tower2'] for k in range(4)]
    def prediction(tau,s,k):return next(x for x in supported if x['tower']==s and x['timescale_reference_days']==tau)['predictions_mmwe'][k]
    def choose(keys):return min(valid_tau,key=lambda tau:(np.mean([abs(prediction(tau,s,k)-truth[s]['development'][k]) for s,k in keys]),tau))
    def choose_calendar(keys):
        gammas=sorted(set(x['parameter'] for x in calendar_pool))
        def p(g,s,k):return next(x for x in calendar_pool if x['tower']==s and x['parameter']==g)['predictions'][k]
        g=min(gammas,key=lambda g:(np.mean([abs(p(g,s,k)-truth[s]['development'][k]) for s,k in keys]),abs(g),g))
        return g,p
    lopo=[]
    for key in cases:
        training=[x for x in cases if x!=key];tau=choose(training);g,cal=choose_calendar(training);s,k=key;y=truth[s]['development'][k]
        lopo.append({'tower':s,'period_index':k,'selected_tau_days_other7':tau,'predicted_mmwe':prediction(tau,s,k),'observed_development_mmwe':float(y),'absolute_error_mmwe':abs(prediction(tau,s,k)-y),'calendar_gamma_other7':g,'calendar_abs_error_mmwe':abs(cal(g,s,k)-y),'locked_air_abs_error_mmwe':abs(fixed_baselines[s]['air_pdd']['predictions_mmwe'][k]-y)})
    loso=[]
    for s in ['Tower1','Tower2']:
        other='Tower2'if s=='Tower1'else'Tower1';keys=[(other,k)for k in range(4)];tau=choose(keys);g,cal=choose_calendar(keys)
        errors=[abs(prediction(tau,s,k)-truth[s]['development'][k])for k in range(4)]
        loso.append({'test_stake_DEVELOPMENT':s,'parameter_stake':other,'selected_tau_days':tau,'MAE_mmwe':float(np.mean(errors)),'errors_mmwe':[float(x)for x in errors],'predictions_mmwe':[prediction(tau,s,k)for k in range(4)],'calendar_gamma_otherstake':g,'calendar_MAE_mmwe':float(np.mean([abs(cal(g,s,k)-truth[s]['development'][k])for k in range(4)])),'locked_air_MAE_mmwe':fixed_baselines[s]['air_pdd']['DEVELOPMENT_MAE_mmwe']})
    report={'scope':'ALLT1/T2targets are alreadyreleased DEVELOPMENT; nofreshvalidation, noT3/laterT1access','family_declared_before_run':{'timescale_reference_days':TIMESCALE_DAYS,'albedo':ALBEDO,'emissivity':EMISSIVITY,'tau_geometry_scaling':'(d/.5)^2','gain':'oneTRAINheatpartitionfraction[0,1]perstake','noatmospheric_turbulence_latent_snowmodel':True},'truth':{s:{'train':v['train'],'development':v['development'].tolist()}for s,v in truth.items()},'all_trials':trials,'locked_old_baselines':fixed_baselines,'retrospective_LOPO':lopo,'retrospective_leave_one_stake_parameter_selection':loso,'pooled_retrospective_LOPO_MAE_mmwe':float(np.mean([x['absolute_error_mmwe']for x in lopo])),'pooled_retrospective_LOPO_calendar_MAE_mmwe':float(np.mean([x['calendar_abs_error_mmwe']for x in lopo])),'pooled_locked_air_MAE_mmwe':float(np.mean([x['locked_air_abs_error_mmwe']for x in lopo])),
    'source_hashes':{str(x):sha(x)for x in [R3/'moisture/DEVELOPMENT_FIT.json',R3/'moisture/frozen_tower2_final/RESULTS.json',R3/'moisture/frozen_tower2_final/PREDICTIONS.json',R3/'data/moisture/PIR22-23_60-kd_ERA5.csv',Path(__file__),Path(__file__).parent/'radiation_reservoir.py']}}
    (output/'DEVELOPMENT_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PooledDEVELOPMENTLOPOreservoir/calendar/air MAE',report['pooled_retrospective_LOPO_MAE_mmwe'],report['pooled_retrospective_LOPO_calendar_MAE_mmwe'],report['pooled_locked_air_MAE_mmwe'])
    print('LOSO',loso)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('output');args=a.parse_args();run(args.output)
