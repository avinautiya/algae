"""Frozen deterministic native-temperature masking; never uses melt outcomes."""
import argparse,bisect,csv,datetime as dt,hashlib,json,math,statistics,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'temperature_digitization'))
from native_temperature import TemperatureSeries,seconds

FROZEN='1a1c5534e1733db621c3912c80dd7fcca001924b2d0ab6c7c905b06ce1545128'

def endpoint(series,time):
    t=seconds(time);j=bisect.bisect_right(series.starts,t)-1
    if j<0 or t>series.ends[j]:raise ValueError('unsupported endpoint')
    a,b,ta,tb=series.segments[j]
    return ta+(tb-ta)*(t-a)/(b-a)

def block_errors(series,start,end):
    truth=series.integrate(start,end)
    if abs(truth['coverage_fraction']-1)>1e-12 or truth['max_uncovered_gap_seconds']>1e-9:
        return None
    left,right=endpoint(series,start),endpoint(series,end);duration=(end-start).total_seconds()
    predictions={'linear_endpoint_bridge':duration*(left+right)/2,
                 'constant_previous_endpoint':duration*left,
                 'constant_next_endpoint':duration*right}
    return [dict(method=m,truth_integral_C_s=truth['integral_C_s'],prediction_integral_C_s=v,
                 signed_integral_error_C_s=v-truth['integral_C_s'],
                 signed_mean_temperature_error_C=(v-truth['integral_C_s'])/duration,
                 true_mean_temperature_C=truth['integral_C_s']/duration,
                 left_temperature_C=left,right_temperature_C=right)
            for m,v in predictions.items()]

def scores(rows):
    e=[r['signed_mean_temperature_error_C'] for r in rows]
    return dict(blocks=len(rows),MAE_C=statistics.mean(abs(x) for x in e),
        MAE_integral_C_s=statistics.mean(abs(r['signed_integral_error_C_s']) for r in rows),
        bias_C=statistics.mean(e),RMSE_C=math.sqrt(statistics.mean(x*x for x in e)),
        median_absolute_error_C=statistics.median(abs(x) for x in e),max_absolute_error_C=max(abs(x) for x in e))

def run(root,out,protocol_path):
    raw=protocol_path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=FROZEN:raise ValueError('frozen protocol changed')
    p=json.loads(raw)
    if hashlib.sha256((root/'SOURCE_MANIFEST.json').read_bytes()).hexdigest()!=p['source_manifest_sha256']:
        raise ValueError('source manifest changed')
    out.mkdir(parents=True,exist_ok=True);allrows=[];counts={};summary={}
    start=dt.datetime.fromisoformat(p['selection_window_start']);end=dt.datetime.fromisoformat(p['selection_window_end_exclusive'])
    for station in p['station_ids']:
        series=TemperatureSeries.from_archive(root,station,p['max_native_gap_seconds'])
        for hours in p['mask_duration_hours']:
            current=start;attempted=0;excluded=[];accepted=0
            while current+dt.timedelta(hours=hours)<=end:
                finish=current+dt.timedelta(hours=hours);attempted+=1
                errors=block_errors(series,current,finish)
                if errors is None:excluded.append(current.isoformat())
                else:
                    accepted+=1
                    for r in errors:r.update(station=f'P{station}',duration_h=hours,start=current.isoformat(),end=finish.isoformat(),id=f'P{station}:{hours}h:{current.isoformat()}')
                    allrows.extend(errors)
                current=finish
            key=f'P{station}_{hours}h'
            counts[key]=dict(candidate_blocks=attempted,retained_blocks=accepted,excluded_incomplete_support_start_times=excluded)
            summary[key]={m:scores([r for r in allrows if r['station']==f'P{station}' and r['duration_h']==hours and r['method']==m])
                          for m in p['methods']} if accepted else {}
    with (out/'MASKED_BLOCK_PREDICTIONS.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
    result=dict(protocol_sha256=FROZEN,source_manifest_sha256=p['source_manifest_sha256'],
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        counts=counts,scores=summary,scope='Controlled measured-floor-temperature reconstruction, not new field model validation',
        notes=p['limitations']+['No duration pooling or independent-block confidence interval.'])
    (out/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('root',type=Path);a.add_argument('--out',type=Path,required=True);a.add_argument('--protocol',type=Path,default=Path(__file__).resolve().parent/'PROTOCOL_FROZEN.json');args=a.parse_args();run(args.root,args.out,args.protocol)
