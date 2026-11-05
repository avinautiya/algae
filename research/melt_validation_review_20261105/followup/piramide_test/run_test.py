"""Frozen, recipient-excluded predictions on camera-read direct ice stakes.

Uses raw Pir AWS fields, never its filled/synthetic tables. The timestamp
interval interpretation is reported in both directions rather than tuned.
"""
from pathlib import Path
from zipfile import ZipFile
from datetime import datetime,timedelta
import json,hashlib
import numpy as np
from decode import table,raw_weather
R=Path(__file__).resolve().parent;P=R.parent
HOUR=3600000

def integrate(time,fields,start,end,label_end):
    left=time-(HOUR if label_end else 0)
    ids=np.where((left<end)&(left+HOUR>start))[0]
    if not len(ids):raise ValueError('no forcing')
    spans=np.maximum(0,np.minimum(left[ids]+HOUR,end)-np.maximum(left[ids],start))/HOUR
    if abs(sum(spans)-(end-start)/HOUR)>1e-8 or (np.diff(left[ids])!=HOUR).any():raise ValueError('forcing clock gap')
    values=np.column_stack([fields[k][ids] for k in ['Ta','SW_in','SW_out']])
    if not np.isfinite(values).all() or (values<=-900).any():raise ValueError('missing measured forcing')
    return values,spans

def main():
    raw=(P/'PIRAMIDE_FREEZE.json').read_bytes();assert hashlib.sha256(raw).hexdigest()==(P/'PIRAMIDE_FREEZE_SHA256.txt').read_text().strip()
    spec=json.loads(raw);models=spec['models'];assert all('Piramide' not in m['donors'] for m in models.values())
    archive=P/'Debris_properties_estimations.zip';source=json.loads((P/'11581073_metadata.json').read_text())['files'][0]
    assert 'md5:'+hashlib.md5(archive.read_bytes()).hexdigest()==source['checksum']
    # Eligibility details fixed before calculating prediction scores. This is a
    # local prospective freeze, not external registration or author escrow.
    eligibility=dict(primary='all complete03:00-to03:00 camera-endpoint days',geometry_m=[.50,.45,.24],forcing='raw Hour.Pir Ta SW_in SW_out with actual timestamps; piecewise-constant hourly means; both start/end label interpretations',secondary='source-code snow-screen stratum: midday radiation ratio <=0.3; all-day primary retained',unit='native distance increment metres ice; ablation column independently satisfies0.917*distance+constant in metres water equivalent despite README mm label',models_fixed_sha256=hashlib.sha256(raw).hexdigest(),limitations=['timestamps and camera quality inherited from source; no new fill','three sites, same glacier campaign','03:00anchor chosen from published camera schedule after decoding table headers/values but before computing scores; not a pristine preregistered endpoint','no independent quantified camera-reading error model','published KP1 original training provenance not fully disjoint from old glacier studies'])
    with (R/'ELIGIBILITY.json').open('x') as f:json.dump(eligibility,f,indent=2)
    predictions=[];refused=[];gates=[]
    with ZipFile(archive) as z:
        wt,wf=raw_weather(z.read('Debris_properties_estimations/Data/Calibration/Meteo_AWS_filled_partition.mat'))
        for site,h in enumerate(eligibility['geometry_m'],1):
            time,cols=table(z.read(f'Debris_properties_estimations/Data/Ablation_Tower{site}.mat'))
            good=np.isfinite(cols['distance'])&np.isfinite(cols['ablation']);identity=cols['ablation'][good]-.917*cols['distance'][good]
            assert np.ptp(identity)<1e-10
            gates.append(dict(site=site,h_m=h,native_records=len(time),finite_stake_records=int(good.sum()),water_equivalent_identity_range=float(np.ptp(identity))))
            indices={int(t):i for i,t in enumerate(time)}
            for i,start in enumerate(time):
                # MATLAB's source datetime is unzoned. Epoch-like milliseconds
                # encode its calendar only; they do not establish physical UTC.
                dt=datetime.utcfromtimestamp(start/1000)
                if dt.hour!=3 or dt.minute!=0:continue
                end=start+24*HOUR;j=indices.get(int(end))
                if j is None or not good[i] or not good[j]:
                    refused.append(dict(site=site,date=dt.isoformat(),reason='missing03h camera endpoints'));continue
                observed=float(cols['distance'][j]-cols['distance'][i])
                for label_end in [False,True]:
                    try:values,spans=integrate(wt,wf,start,end,label_end)
                    except ValueError as exc:
                        refused.append(dict(site=site,date=dt.isoformat(),label_end=label_end,reason=str(exc)));continue
                    temp=float(spans@np.maximum(values[:,0],0));solar=float(spans@np.maximum(values[:,1]-values[:,2],0));meanT=float(spans@values[:,0]/24)
                    d=dict(constant=24,pdd=temp,shortwave=solar,pdd_over_h=temp/h,shortwave_over_h=solar/h)
                    # Match the source-code midday albedo criterion. This is
                    # still a proxy, not independent confirmation of no snow.
                    noon_start=start+(11-3)*HOUR
                    mids=np.where((wt>=noon_start)&(wt<noon_start+3*HOUR))[0]
                    ratios=wf['SW_out'][mids]/wf['SW_in'][mids]
                    ratio=float(np.mean(ratios)) if len(mids)==3 and np.isfinite(ratios).all() and ((ratios>=0)&(ratios<=1)).all() else float('nan')
                    issued={}
                    for name,m in models.items():issued[name]=float(np.array([d[f] for f in m['features']])/m['feature_scale']@np.array(m['beta'])*m['target_scale'])
                    # Winter-Billington2020 Table4 fixed effects, zero random
                    # effect for an unseen glacier: mmwe/day/K -> metres ice.
                    issued['published_KP1_2020_fixed_effects']=10**(.62-1.46*h)*max(meanT,0)/917
                    for name,p in issued.items():predictions.append(dict(method=name,site=site,date=dt.isoformat(),end=datetime.utcfromtimestamp(end/1000).isoformat(),label_end=label_end,h_m=h,observed_m_ice=observed,predicted_m_ice=p,error_m=p-observed,albedo_proxy=ratio if np.isfinite(ratio) else None,source_snow_screen=bool(np.isfinite(ratio) and ratio<=.3),features=d,mean_daily_T=meanT))
    scores=[]
    for label_end in [False,True]:
        for stratum in ['all_complete_days','albedo_proxy_le_0.3']:
            for name in [*models,'published_KP1_2020_fixed_effects']:
                rows=[r for r in predictions if r['method']==name and r['label_end']==label_end and (stratum=='all_complete_days' or r['source_snow_screen'])]
                per_site={}
                for site in [1,2,3]:
                    sub=[r for r in rows if r['site']==site]
                    if not sub:continue
                    errors=np.array([r['error_m'] for r in sub])
                    per_site[str(site)]=dict(n=len(sub),MAE_m=float(np.mean(abs(errors))),bias_m=float(np.mean(errors)),observed_sum_m=sum(r['observed_m_ice'] for r in sub),predicted_sum_m=sum(r['predicted_m_ice'] for r in sub))
                if per_site:scores.append(dict(method=name,label_end=label_end,stratum=stratum,sites=per_site,equal_site_MAE_m=float(np.mean([s['MAE_m'] for s in per_site.values()]))))
    assert predictions
    result=dict(model_freeze_sha256=hashlib.sha256(raw).hexdigest(),archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),eligibility_sha256=hashlib.sha256((R/'ELIGIBILITY.json').read_bytes()).hexdigest(),gates=gates,refused=refused,scores=scores,rows=predictions,scope='direct camera-stake sub-debris ice ablation; coefficients frozen, eligible-day rules determined during observation audit; observed meteorology and uncertain interval labels',limitations=eligibility['limitations'])
    with (R/'RESULTS.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps([dict(method=s['method'],stratum=s['stratum'],label_end=s['label_end'],MAE_mm=s['equal_site_MAE_m']*1000,n_by_site={k:v['n'] for k,v in s['sites'].items()}) for s in scores],indent=2))

if __name__=='__main__':main()
