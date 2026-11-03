"""Target-free forcing diagnostics and melt-to-energy dimensional accounting.

No fit, field target, or melt prediction. Radiation is measured AWS surface flux;
it is NOT the conductive flux delivered to ice beneath debris.
"""
import csv,datetime as dt,hashlib,io,json,math,zipfile
from pathlib import Path
SIGMA=5.670374419e-8
LATENT=334000.
RHOW=1000.

def equivalent_flux(delta_mwe,seconds):
    if not all(type(x) in (int,float) and math.isfinite(x) for x in (delta_mwe,seconds)) or seconds<=0:
        raise ValueError('finite melt difference and positive duration required')
    return RHOW*LATENT*delta_mwe/seconds

def radiometric_temperature(lwout,lwin,eps):
    if not all(math.isfinite(x) for x in (lwout,lwin,eps)) or not 0<eps<=1:raise ValueError('invalid radiometry')
    numerator=lwout-(1-eps)*lwin
    if numerator<=0:raise ValueError('nonphysical radiometry')
    return (numerator/(eps*SIGMA))**.25-273.15

def process_window(rows,start,end,eps=.93):
    # Date-label support only, NOT a claim about unknown subdaily stake visit time.
    rr=[r for r in rows if start<=r['time']<end]
    wanted=int((end-start).total_seconds()/3600)
    if len(rr)!=wanted or any((b['time']-a['time']).total_seconds()!=3600 for a,b in zip(rr,rr[1:])):
        raise ValueError('complete exact hourly forcing support required')
    for r in rr:
        if not all(math.isfinite(r[k]) and r[k]!=-9999 for k in ('T','FF','RH','SWIN','SWOUT','LWIN','LWOUT','PP','SNOW')):
            raise ValueError('invalid forcing; no silent fill')
    mean=lambda k:sum(k(r) for r in rr)/len(rr)
    return dict(start_label=start.isoformat(),end_label_exclusive=end.isoformat(),hours=len(rr),
      mean_air_temperature_C=mean(lambda r:r['T']),mean_positive_air_temperature_C=mean(lambda r:max(r['T'],0.)),mean_wind_m_s=mean(lambda r:r['FF']),
      mean_net_shortwave_W_m2=mean(lambda r:r['SWIN']-r['SWOUT']),
      mean_net_longwave_W_m2=mean(lambda r:r['LWIN']-r['LWOUT']),
      mean_net_radiation_W_m2=mean(lambda r:r['SWIN']-r['SWOUT']+r['LWIN']-r['LWOUT']),
      precipitation_label_hours=sum(r['PP']>0 for r in rr),
      positive_snow_label_hours=sum(r['SNOW']>0 for r in rr),
      dry_label_hours_with_RH_below_50_percent=sum(r['PP']==0 and r['RH']<50 for r in rr),
      mean_radiometric_surface_temperature_assumed_eps_C=mean(lambda r:radiometric_temperature(r['LWOUT'],r['LWIN'],eps)),
      air_warmer_than_surface_hours_assumed_eps=sum(r['T']>radiometric_temperature(r['LWOUT'],r['LWIN'],eps) for r in rr),
      assumed_emissivity=eps,independent_surface_temperature_measurement=False,
      actual_evaporation_measured=False,ice_flux_inferred_from_surface_radiation=False)

def main(source,out):
    raw=source.read_bytes()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        name='PIR_20142015/PIR_20142015_60.csv';text=z.read(name)
    rows=[]
    for x in csv.DictReader(io.StringIO(text.decode())):
        hour=int(float(x['HOUR'])); minute=int(float(x['MINUTE']))
        if not 0<=hour<=24 or not 0<=minute<60 or (hour==24 and minute!=0):raise ValueError('invalid source clock')
        t=dt.datetime(*(int(float(x[k])) for k in ['YEAR','MONTH','DAY']))+dt.timedelta(hours=hour,minutes=minute)
        rows.append(dict(time=t,**{k:float(x[k]) for k in ['T','RH','FF','SWIN','SWOUT','LWIN','LWOUT','PP','SNOW']}))
    windows=[('2014-12-04','2015-02-10'),('2015-02-10','2015-04-08')]
    result=dict(source_archive_sha256=hashlib.sha256(raw).hexdigest(),source_member_sha256=hashlib.sha256(text).hexdigest(),
      output_scope='Source-label forcing diagnostics only; no target or model score read.',
      windows=[process_window(rows,dt.datetime.fromisoformat(a),dt.datetime.fromisoformat(b)) for a,b in windows],
      generic_dimensional_reference_not_observed_error=dict(delta_ice_ablation_mwe=.1,duration_days=30,mean_heat_difference_W_m2=equivalent_flux(.1,30*86400)),
      causal_attribution_established=False,novel_melt_model_established=False)
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('out',type=Path);a=p.parse_args();main(a.source,a.out)
