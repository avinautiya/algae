"""Positive empirical melt indices; no physical energy-conservation assertion.

The primary is invariant to a constant additive moisture-sensor offset after
local calibration. Other centered functions can share that offset invariance.
No causal moisture law, historical novelty, or universal accuracy is implied.
"""
import csv,datetime as dt,hashlib,json,math,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from paths import require_piramide
ROOT=Path(__file__).resolve().parents[1]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def positive_area(a,b,seconds):
    if min(a,b)>=0:return seconds*(a+b)/2
    if max(a,b)<=0:return 0.
    return seconds*max(a,b)**2/(2*abs(b-a))

def air_exposure(labels,source):
    with Path(source).open()as f:rows=list(csv.DictReader(f))
    origin=labels[0];t=np.array([(dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL'])-origin).total_seconds()for r in rows]);v=np.array([float(r['TA_3'])for r in rows])
    if not np.isfinite(v).all()or np.any(np.diff(t)<=0):raise ValueError('invalid native air stream')
    boundaries=np.array([(x-origin).total_seconds()for x in labels]);out=[];largest=0.
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        if a<t[0]or b>t[-1]:raise ValueError('no air extrapolation permitted')
        i=max(int(np.searchsorted(t,a,side='right'))-1,0);total=0.
        while i<len(t)-1 and t[i]<b:
            left,right=max(a,t[i]),min(b,t[i+1]);gap=t[i+1]-t[i]
            if right>left:
                if gap>64800:raise ValueError('native air gap exceeds declared18h limit')
                largest=max(largest,gap)
                vl=v[i]+(v[i+1]-v[i])*(left-t[i])/gap;vr=v[i]+(v[i+1]-v[i])*(right-t[i])/gap
                total+=positive_area(vl,vr,right-left)
            i+=1
        out.append(total)
    return np.asarray(out),largest

def load_site(tower):
    if tower not in ['Tower1','Tower2']:raise ValueError('declared tower required')
    folder=ROOT/'data/moisture';required=json.loads((folder/'REQUIRED_FREEZE_FILES.json').read_text())['files']
    for name in [tower+'_CANONICAL_FEATURES.csv',tower+'_NATIVE_AIR_TA3.csv','CANONICAL_FEATURE_MANIFEST.json']:
        if sha(folder/name)!=required[name]:raise ValueError('canonical feature source changed')
    with (folder/(tower+'_CANONICAL_FEATURES.csv')).open()as f:rows=list(csv.DictReader(f))
    labels=[dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL'])for r in rows]
    if any((b-a).total_seconds()!=1800 for a,b in zip(labels,labels[1:])):raise ValueError('complete native halfhour feature support required')
    theta=np.array([float(r['SM3_raw_m3_m3'])for r in rows[:-1]])
    if not np.isfinite(theta).all()or np.min(theta+.031)<=0:raise ValueError('shifted-ratio input unsupported; no clipping or sample deletion')
    air,maxgap=air_exposure(labels,folder/(tower+'_NATIVE_AIR_TA3.csv'))
    previous=require_piramide();temperature_path=previous/'data/decoded'/('debtemp_'+tower+'.csv')
    manifest=json.loads((previous/'data/decoded/DECODE_VERIFICATION.json').read_text())
    expected=next(x['sha256']for x in manifest['checks']if x['family']=='debtemp'and x['tower']==tower)
    if sha(temperature_path)!=expected:raise ValueError('source temperature table changed')
    with temperature_path.open()as f:temp={r['SOURCE_LOCAL_LABEL']:r for r in csv.DictReader(f)}
    trows=[temp[x.isoformat()]for x in labels[:-1]]
    upper=np.array([float(r['DT1'])for r in trows]);spacing=.125 if tower=='Tower1'else .11
    gradient=np.array([(float(r['DT4'])-float(r['DT5']))/spacing for r in trows])
    if not np.isfinite(upper).all()or not np.isfinite(gradient).all():raise ValueError('incomplete profile feature support')
    epochs=json.loads((folder/'CANONICAL_FEATURE_MANIFEST.json').read_text())['sites'][tower]['frozen_site_specific_anchors']
    endpoints=[labels.index(dt.datetime.fromisoformat(x))for x in epochs]
    slices=[slice(a,b)for a,b in zip(endpoints[:-1],endpoints[1:])]
    reference=float(np.mean(theta[slices[0]]));days=np.array([((x-labels[0]).total_seconds()+900)/86400 for x in labels[:-1]])-7.
    return dict(tower=tower,labels=labels,theta=theta,theta_reference=reference,air=air,
        buried=np.maximum(upper,0)*1800,gradient=np.maximum(gradient,0)*1800,
        seconds=np.full(len(theta),1800.),days=days,slices=slices,epochs=epochs,
        source_max_air_gap_seconds=float(maxgap),soil_subzero_reference_rows=sum(r['soilT_reference_below0']=='True'for r in rows[:-1]),
        metadata='Empirical indices: positive gradient feature is not a signed physical heat flux; no energy or manufacturer-precision certification.')

def dose(features,method,parameter=0.):
    f=features;z=f['theta']-f['theta_reference'];p=float(parameter)
    if not math.isfinite(p):raise ValueError('finite parameter required')
    if method=='anomaly_exp_air':values=f['air']*np.exp(p*z)
    elif method=='air_pdd':values=f['air']
    elif method=='buried_pdd':values=f['buried']
    elif method=='gradient_index':values=f['gradient']
    elif method=='constant_rate':values=f['seconds']
    elif method=='shifted_ratio_air':values=f['air']*((f['theta']+.031)/(f['theta_reference']+.031))**p
    elif method=='centered_softplus_air':values=f['air']*np.logaddexp(0.,p*z)
    elif method=='calendar_exp_air':values=f['air']*np.exp(p*f['days'])
    elif method=='moisture_only_exp':values=f['seconds']*np.exp(p*z)
    else:raise ValueError('unknown predeclared index')
    if not np.isfinite(values).all()or np.min(values)<0:raise ValueError('nonfinite or negative dose; no clipping')
    return np.array([math.fsum(values[s])for s in f['slices']])

def prediction_ratios(features,method,parameter):
    h=dose(features,method,parameter)
    if h[0]<=0:raise ValueError('positive train exposure required')
    return h[1:]/h[0]

def development_fit(features,truth_train,truth_development,method,grid):
    y=np.asarray(truth_development,float)
    rows=[]
    for p in grid:
        pred=truth_train*prediction_ratios(features,method,p);loss=float(np.mean(np.abs(pred-y)))
        rows.append(dict(parameter=float(p),development_MAE_mmwe=loss,predictions_mmwe=pred.tolist()))
    winner=min(rows,key=lambda r:(r['development_MAE_mmwe'],abs(r['parameter']),r['parameter']))
    loo=[]
    for hold in range(len(y)):
        selected=min(rows,key=lambda r:(float(np.mean(np.delete(np.abs(np.array(r['predictions_mmwe'])-y),hold))),abs(r['parameter']),r['parameter']))
        loo.append(dict(held_period=hold,parameter=selected['parameter'],absolute_error_mmwe=abs(selected['predictions_mmwe'][hold]-float(y[hold]))))
    return dict(method=method,selected=winner,grid=rows,selected_at_search_boundary=winner['parameter']in[float(grid[0]),float(grid[-1])],
        leave_one_development_period_out_diagnostic=loo,new_validation_claim=False)

def model_grids():
    return dict(anomaly_exp_air=np.linspace(-40,40,801),shifted_ratio_air=np.linspace(-4,4,801),
        centered_softplus_air=np.linspace(-80,80,801),calendar_exp_air=np.linspace(-.04,.04,801),moisture_only_exp=np.linspace(-40,40,801))

def develop(output):
    if output.exists():raise ValueError('new development output required')
    previous=require_piramide();p=json.loads((previous/'model/frozen_piramide_v3/melt_comparison/PREDICTIONS.json').read_text())
    first=json.loads((previous/'model/frozen_piramide_v3/melt_comparison/RESULTS.json').read_text());batch=json.loads((previous/'model/frozen_piramide_future_batch_final/melt_comparison/RESULTS.json').read_text())
    name=first['primary_candidate'];y0=p['calibration']['observed_ice_ablation_mmwe'];y=[first['model_metrics'][name]['observed_total_mmwe']]+[r['observed_mmwe']for r in batch['model_metrics'][name]['period_errors']]
    f=load_site('Tower1');fits={m:development_fit(f,y0,y,m,g)for m,g in model_grids().items()}
    result=dict(scope='All four Tower1 periods are now development data for this NEW method; prior failed results preserved',
        train_mmwe=y0,development_mmwe=y,epochs=f['epochs'],fits=fits,primary='anomaly_exp_air',
        Tower2_camera_values_opened=False,new_field_validation=False,
        source_prior_predictions_sha256=sha(previous/'model/frozen_piramide_v3/melt_comparison/PREDICTIONS.json'),
        source_prior_first_result_sha256=sha(previous/'model/frozen_piramide_v3/melt_comparison/RESULTS.json'),
        source_prior_batch_result_sha256=sha(previous/'model/frozen_piramide_future_batch_final/melt_comparison/RESULTS.json'))
    with output.open('x')as o:json.dump(result,o,indent=2);o.write('\n')
    return {m:dict(parameter=x['selected']['parameter'],DEVELOPMENT_MAE=x['selected']['development_MAE_mmwe'])for m,x in fits.items()}

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('output',type=Path);a=p.parse_args();print(json.dumps(develop(a.output),indent=2))
