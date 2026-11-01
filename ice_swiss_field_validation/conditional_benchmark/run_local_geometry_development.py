"""V5 target-exposed local-hydraulics development; no test-based tuning."""
import argparse
import csv
import datetime as dt
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent/'temperature_digitization'))
from native_temperature import TemperatureSeries
from standard_closures import nusselt

MODELS = ['constant_rate_fitted', 'constant_h_fitted', 'db04_unfitted',
          'db03_unfitted', 'gnielinski_smooth_unfitted', 'db04_scale_fitted']
K, PR, RHO, LATENT = .56, 13.5, 917., 334000.

def timestamp(s):
    return dt.datetime(*[int(x) for x in re.findall(r'\d+', s)[:6]])

def filehash(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

class AuthorDailySeries:
    """Explicit scenario, not a native observed full-day temperature integral."""
    def __init__(self, source, station):
        with (source/f'products/temp{station}_daily.csv').open() as f:
            self.days={timestamp(r['Date_Time']).date():float(r['Temperature_daily']) for r in csv.DictReader(f)}
    def integrate_positive(self, start, end, extra_offset_C=0):
        cursor=start;integral=covered=0.;largest=0.
        while cursor<end:
            stop=min(end,cursor.replace(hour=0,minute=0,second=0)+dt.timedelta(days=1))
            seconds=(stop-cursor).total_seconds();value=self.days.get(cursor.date(),float('nan'))
            if math.isfinite(value):
                integral+=max(0.,value+extra_offset_C)*seconds;covered+=seconds
            else:largest=max(largest,seconds)
            cursor=stop
        return dict(integral_C_s=integral,covered_seconds=covered,requested_seconds=(end-start).total_seconds(),max_uncovered_gap_seconds=largest)

def targets(source):
    allrows=[]
    for station in [3,5]:
        path=source/f'data/Channel_elevation/P{station}_channel_elevation.csv'
        rows=list(csv.reader(path.open(), delimiter=';'))
        epochs=[(timestamp(r[0]),float(r[1]),float(r[2])) for r in rows[1:] if r]
        split=max(t for t,z,e in epochs if t.date()<=dt.date(2019,7,25))
        for j,(a,b) in enumerate(zip(epochs,epochs[1:])):
            role='train' if b[0]<=split else 'test'
            if role=='test' and a[0]<split:raise AssertionError('split crosses interval')
            allrows.append(dict(id=f'P{station}:{j}',station=f'P{station}',role=role,
                                start=a[0].isoformat(),end=b[0].isoformat(),
                                observed_lowering_m=a[1]-b[1],epoch0_std_m=a[2],epoch1_std_m=b[2],
                                summed_epoch_error_scenario_m=a[2]+b[2]))
    return allrows

def fit_scale(rows, column):
    denominator=math.fsum(r[column]**2 for r in rows)
    if denominator<=0:raise ValueError('training basis has zero information')
    return max(0.,math.fsum(r[column]*r['observed_lowering_m'] for r in rows)/denominator)

def metrics(rows, name):
    e=[r['predictions_m'][name]-r['observed_lowering_m'] for r in rows]
    if not e:return dict(n=0)
    return dict(n=len(e),MAE_m=math.fsum(abs(x) for x in e)/len(e),
                RMSE_m=math.sqrt(math.fsum(x*x for x in e)/len(e)),bias_m=math.fsum(e)/len(e),
                observed_total_m=math.fsum(r['observed_lowering_m'] for r in rows),
                predicted_total_m=math.fsum(r['predictions_m'][name] for r in rows))

def run(source, output, forcing='native', temperature_factory=None):
    freeze=json.loads((ROOT/'FREEZE_V2.json').read_text())
    if freeze['protocol_sha256']!=filehash(ROOT/'PROTOCOL_V2_ORIGINAL_DATA.md'):
        raise ValueError('protocol changed after freeze')
    if freeze['source_manifest_sha256']!=filehash(source/'SOURCE_MANIFEST.json'):
        raise ValueError('different original-data manifest')
    if forcing=='author_daily':
        secondary=json.loads((ROOT/'FREEZE_DAILY_SCENARIO.json').read_text())
        if secondary['secondary_protocol_sha256']!=filehash(ROOT/'SECONDARY_DAILY_SCENARIO.md'):
            raise ValueError('secondary protocol changed')
        freeze=dict(freeze,secondary=secondary)
        series={f'P{s}':AuthorDailySeries(source,s) for s in [3,5]}
    elif temperature_factory is not None:
        series={f'P{s}':temperature_factory(source,s) for s in [3,5]}
    else:
        series={f'P{s}':TemperatureSeries.from_archive(source,s,max_gap_seconds=60) for s in [3,5]}
    with (source/'products/Reynold_daily.csv').open() as f:
        reynolds={timestamp(r['Date_Time']).date():float(r['Reynold']) for r in csv.DictReader(f)}
    rawtargets=targets(source)
    def daily(filename,column):
        with (source/'products'/filename).open() as stream:
            return {timestamp(row['Date_Time']).date():float(row[column]) for row in csv.DictReader(stream)}
    discharge=daily('discharge_channel_full_daily.csv','Discharge_daily')
    local_depth={'P3':daily('stage3_full_daily.csv','Water_height'), 'P5':daily('stage5_daily.csv','Water_height')}
    @lru_cache(maxsize=None)
    def temperature_integral(station,start,end,shift):
        return series[station].integrate_positive(start,end,extra_offset_C=shift)

    @lru_cache(maxsize=None)
    def basis(station,start,end,clock_shift,temperature_shift,width):
        a=timestamp(start)+dt.timedelta(hours=clock_shift)
        b=timestamp(end)+dt.timedelta(hours=clock_shift)
        duration=(b-a).total_seconds();cells=[];cursor=a;covered=0.;maxgap=0.
        while cursor<b:
            nextday=cursor.replace(hour=0,minute=0,second=0)+dt.timedelta(days=1)
            stop=min(b,nextday)
            flow=discharge.get(cursor.date(),float('nan'));depth=local_depth[station].get(cursor.date(),float('nan'))
            if not math.isfinite(flow) or not math.isfinite(depth) or flow<=0 or depth<=0:
                return {'excluded':'missing or nonpositive imposed local daily depth/discharge'}
            diameter=4*width*depth/(width+2*depth)
            velocity=flow/(width*depth);reyn=velocity*diameter/1.8e-6
            try:
                nu={n:nusselt(reyn,PR,m) for n,m in [('db04','db_published_04'),('db03','db_cooling_03'),('gn','gnielinski_smooth')]}
            except ValueError as exc:return {'excluded':str(exc)}
            integral=temperature_integral(station,cursor,stop,temperature_shift)
            covered+=integral['covered_seconds'];maxgap=max(maxgap,integral['max_uncovered_gap_seconds'])
            cells.append(dict(start=cursor.isoformat(),end=stop.isoformat(),Re=reyn,Nu=nu,local_depth_m=depth,discharge_m3_s=flow,velocity_m_s=velocity,hydraulic_diameter_m=diameter,channel_width_m=width,
                              integral_C_s=integral['integral_C_s']))
            cursor=stop
        if duration-covered>1e-5:
            return {'excluded':'incomplete selected temperature forcing support',
                    'coverage':dict(requested_seconds=duration,covered_seconds=covered,
                                    coverage_fraction=covered/duration,max_daily_uncovered_gap_seconds=maxgap)}
        return dict(duration_s=duration,temperature_integral_C_s=math.fsum(c['integral_C_s'] for c in cells),cells=cells)

    scenarios={}
    for clock_shift in [-2,0,2]:
        for tshift in [-.05,0.,.05]:
            for diameter in [1.5,2.0,2.5]:
                name=f'clock{clock_shift:+d}h_T{tshift:+.2f}C_W{diameter:.2f}m'
                accepted=[];excluded=[]
                for target in rawtargets:
                    x=basis(target['station'],target['start'],target['end'],clock_shift,tshift,diameter)
                    if 'excluded' in x:
                        excluded.append(dict(target,reason=x));continue
                    record=dict(target,duration_s=x['duration_s'],temperature_integral_C_s=x['temperature_integral_C_s'],hydraulic_cells=x['cells'])
                    for n in ['db04','db03','gn']:
                        record[n+'_basis_m']=K/(RHO*LATENT)*math.fsum(c['Nu'][n]/c['hydraulic_diameter_m']*c['integral_C_s'] for c in x['cells'])
                    record['constant_h_basis_m']=x['temperature_integral_C_s']/(RHO*LATENT)
                    accepted.append(record)
                train=[x for x in accepted if x['role']=='train'];test=[x for x in accepted if x['role']=='test']
                if not train or not test:
                    scenarios[name]=dict(status='insufficient common-support training or test',excluded_intervals=excluded);continue
                rates=fit_scale(train,'duration_s');h=fit_scale(train,'constant_h_basis_m');scale=fit_scale(train,'db04_basis_m')
                for row in accepted:
                    row['predictions_m']=dict(zip(MODELS,[rates*row['duration_s'],h*row['constant_h_basis_m'],
                        row['db04_basis_m'],row['db03_basis_m'],row['gn_basis_m'],scale*row['db04_basis_m']]))
                scores={site:{model:metrics([r for r in test if site=='pooled' or r['station']==site],model) for model in MODELS} for site in ['P3','P5','pooled']}
                scenarios[name]=dict(status='completed',clock_shift_h=clock_shift,temperature_shift_C=tshift,
                    channel_width_m=diameter,calibration=dict(training_ids=[r['id'] for r in train],
                    constant_rate_m_s=rates,constant_h_W_m2_K=h,db04_scale=scale),
                    test_metrics=scores,accepted_intervals=accepted,excluded_intervals=excluded)
    result=dict(scope='conditional thermal-component retrospective temporal test; one Swiss stream event',
        freeze=freeze,code_sha256=filehash(Path(__file__)),input_targets=rawtargets,
        reconstruction=('native piecewise-linear clipped nonnegative Tw' if forcing=='native' else 'explicit scenario: author sample-weighted daily Tw imposed constant throughout calendar day; not fully observed heat integral')+'; author daily Re held constant within calendar day',
        properties=dict(k_W_m_K=K,Pr=PR,rho_i_kg_m3=RHO,latent_J_kg=LATENT,nu_m2_s=1.8e-6),
        caveats=['floor temperature is not independently observed bulk temperature',
                 'daily local depths and conserved downstream discharge include reconstructed hydraulic forcing',
                 'clock and temperature/diameter scenarios are engineering scenarios, not confidence limits',
                 'correlated interval endpoints and one stream prevent general consistency claims'],scenarios=scenarios)
    output.mkdir(parents=True,exist_ok=True)
    (output/'CONDITIONAL_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    nominal=scenarios['clock+0h_T+0.00C_W2.00m']
    print(json.dumps(dict(status=nominal['status'],calibration=nominal.get('calibration'),test_metrics=nominal.get('test_metrics'),excluded=len(nominal.get('excluded_intervals',[]))),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--forcing',choices=['native','author_daily'],default='native');a=p.parse_args();run(a.source,a.out,a.forcing)
