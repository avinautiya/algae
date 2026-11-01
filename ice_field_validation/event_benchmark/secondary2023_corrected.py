"""Secondary year field hindcast; no melt parameter refitting or winner selection."""
import hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd
import run as base

ROOT=Path(__file__).resolve().parent
MODELS=['training_pooled_constant_normal_rate','training_positive_degree_day',
        'uncalibrated_neutral_energy_balance','training_calibrated_neutral_energy_balance']

def normalize_pair(target,amazon,pairs):
    merged=target.merge(amazon,on='Date_Time',suffixes=('_EB','_AMAZON'),validate='one_to_one')
    normalized={}
    for y,x in pairs.items():
        xcol=x+'_AMAZON' if x in target.columns else x
        ycol=y+'_EB' if y in amazon.columns else y
        normalized[y+'__x']=pd.to_numeric(merged[xcol],errors='coerce')
        normalized[y+'__y']=pd.to_numeric(merged[ycol],errors='coerce')
    return pd.DataFrame(normalized)

def evaluate(rows,w,cal,sites):
    support=base.Clock(w['centers'],np.ones(len(w['centers'])));events=[];excluded=[];clocks={}
    for r in rows:
        if r['site'] not in sites:continue
        ident=r['site']+':'+r['row0']+'-'+r['row1']
        reason=None
        if abs(float(r['dx_horizontal_m'])-.1)>1e-9:reason='not exact10cm marker'
        elif r['t0']<support.edges[0] or r['t1']>support.edges[-1]:reason='outside forcing support'
        if reason:excluded.append(dict(id=ident,reason=reason));continue
        if r['site'] not in clocks:
            sw=base.face(w,float(r['slope_deg']),float(r['aspect_deg']))
            rate_arrays=[np.full(len(sw),cal['constant_normal_rate_m_s']),
                         cal['degree_factor_m_per_Kday']*np.maximum(w['Ta'],0)/86400,
                         base.rate(w,sw,(.2,.002)),base.rate(w,sw,(cal['albedo'],cal['neutral_exchange']))]
            clocks[r['site']]={name:base.Clock(w['centers'],rate) for name,rate in zip(MODELS,rate_arrays)}
        target=.1*math.sin(math.radians(float(r['slope_deg'])))
        predicted={name:clock.arrival(r['t0'],target) for name,clock in clocks[r['site']].items()}
        events.append(dict(id=ident,site=r['site'],observed_duration_h=(r['t1']-r['t0'])/3600,
                           predicted_duration_h={k:(v-r['t0'])/3600 if v is not None else None for k,v in predicted.items()}))
    metrics={}
    for site in sites+['pooled']:
        allrows=[r for r in events if site=='pooled' or r['site']==site]
        common=[r for r in allrows if all(v is not None for v in r['predicted_duration_h'].values())]
        if not common:metrics[site]=dict(eligible=len(allrows),paired_complete=0);continue
        metrics[site]=dict(eligible=len(allrows),paired_complete=len(common),models={});baseline=None
        obs=np.array([r['observed_duration_h'] for r in common])
        for name in MODELS:
            errors=np.array([r['predicted_duration_h'][name] for r in common])-obs;mae=float(np.mean(abs(errors)))
            if baseline is None:baseline=mae
            metrics[site]['models'][name]=dict(MAE_h=mae,RMSE_h=float(np.sqrt(np.mean(errors**2))),bias_h=float(errors.mean()),
                MAPE=float(np.mean(abs(errors)/obs)),fraction_within_2h=float(np.mean(abs(errors)<=2)),
                fraction_within20percent=float(np.mean(abs(errors)/obs<=.2)),MAE_skill_vs_training_constant=1-mae/baseline,
                censored=sum(r['predicted_duration_h'][name] is None for r in allrows))
    return dict(metrics=metrics,events=events,exclusions=excluded)

def transfer_weather(protocol):
    amazon=pd.read_csv(base.RAW/'AWS/Amazon_Creek.csv');pairs=protocol['weather_pairs']
    joined=[]
    for year in (2021,2022):
        target=pd.read_csv(base.RAW/'AWS'/f'Debris_EB_{year}.csv')
        joined.append(normalize_pair(target,amazon,pairs))
    train=pd.concat(joined,ignore_index=True)
    begin,end=protocol['forcing_window']
    year2023=amazon[(amazon.Date_Time>=begin)&(amazon.Date_Time<end)].copy()
    derived=pd.DataFrame(dict(Date_Time=year2023.Date_Time));coefficients={};qc={}
    for target,source in pairs.items():
        x=train[target+'__x'].to_numpy();y=train[target+'__y'].to_numpy()
        valid=np.isfinite(x)&np.isfinite(y)
        a,b=np.linalg.lstsq(np.c_[np.ones(valid.sum()),x[valid]],y[valid],rcond=None)[0]
        raw=a+b*pd.to_numeric(year2023[source],errors='coerce').to_numpy()
        if not np.isfinite(raw).all():raise ValueError('missing secondary forcing '+source)
        clipped=np.clip(raw,0,100) if target=='RH_2m' else (np.maximum(raw,0) if target in ('Wspd','SWin','LWin') else raw)
        derived[target]=clipped;qc[target]=int(np.count_nonzero(raw!=clipped))
        coefficients[target]=dict(source=source,intercept=float(a),slope=float(b),paired_rows=int(valid.sum()),
                                  calibration_weather_RMSE=float(np.sqrt(np.mean((a+b*x[valid]-y[valid])**2))))
    return derived,coefficients,qc

def main():
    protocol=json.loads((ROOT/'protocol2023.json').read_text());cal=json.loads((ROOT/'CALIBRATION_2021.json').read_text())
    frame,coeff,qc=transfer_weather(protocol)
    (ROOT/'WEATHER_TRANSFER_FROZEN.json').write_text(json.dumps(dict(coefficients=coeff,qc=qc,
        protocol_sha256=hashlib.sha256((ROOT/'protocol2023.json').read_bytes()).hexdigest()),indent=2)+'\n')
    # Test targets are read only AFTER weather transfer is saved.
    rows=base.intervals(2023);out={}
    for name,factor in (('base',1.),('SW_minus15percent',.85),('SW_plus15percent',1.15)):
        f=frame.copy();f['SWin']*=factor
        out[name]=evaluate(rows,base.weather(2023,f),cal,protocol['test_sites'])
    result=dict(evaluation='secondary retrospective field hindcast; uncertain offglacier meteorological transfer',
                original_melt_calibration=cal,weather_transfer=coeff,weather_qc_clips=qc,scenarios=out,
                limitations=protocol,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'FIELD_RESULTS_2023.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(transfer=coeff,metrics=out['base']['metrics'],sensitivity_pooled={k:v['metrics']['pooled'] for k,v in out.items()}),indent=2))
if __name__=='__main__':main()
