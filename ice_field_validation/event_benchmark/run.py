"""Held-out-year event-time hindcast. Not a coupled hydraulic-model validation."""
import csv,datetime,hashlib,json,math,os
from pathlib import Path
import numpy as np
import pandas as pd
import pvlib
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parent
DATA=Path(os.environ.get('ALGAE_FIELD_DATA_ROOT',str(ROOT.parent/'kennicott')))
RAW=DATA/'Petersen_etal_2025_Kennicott_Glacier_Data'
RHO=917.;LATENT=334000.;TM=273.15

class Clock:
    def __init__(self,centers,rate):
        self.centers=np.asarray(centers,float);self.rate=np.asarray(rate,float)
        if len(self.centers)!=len(self.rate) or len(self.rate)<2 or not np.isfinite(self.rate).all() or np.any(self.rate<0):raise ValueError('bad nonnegative rate series')
        if not np.allclose(np.diff(self.centers),3600.,rtol=0,atol=.01):raise ValueError('missing/irregular hourly forcing')
        self.edges=np.r_[self.centers-1800,self.centers[-1]+1800]
        self.cum=np.r_[0.,np.cumsum(self.rate*3600)]
    def at(self,t):
        if t<self.edges[0] or t>self.edges[-1]:raise ValueError('outside forcing support')
        j=np.searchsorted(self.edges,t,side='right')-1
        return self.cum[-1] if j==len(self.rate) else self.cum[j]+self.rate[j]*(t-self.edges[j])
    def integral(self,a,b):return self.at(b)-self.at(a)
    def arrival(self,start,target):
        if target<=0:raise ValueError('positive target required')
        goal=self.at(start)+target
        if goal>self.cum[-1]:return None
        j=int(np.searchsorted(self.cum,goal,side='left'))
        if j==0:return float(self.edges[0])
        if self.rate[j-1]<=0:raise RuntimeError('inconsistent inverse cumulative flux')
        return float(self.edges[j-1]+(goal-self.cum[j-1])/self.rate[j-1])

def stamp(s):return pd.Timestamp(s,tz='UTC').timestamp()

def intervals(year):
    geom={}
    base=RAW/'Ice_Cliff_Backwasting'
    with (base/'ice_cliff_slope.csv').open(encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            if row.get(str(year)):geom[row['Stake']]=dict(slope_deg=row[str(year)])
    with (base/'ice_cliff_aspect.csv').open(encoding='utf-8-sig') as stream:
        for row in csv.DictReader(stream):
            if row.get(str(year)) and row['Stake'] in geom:geom[row['Stake']]['aspect_deg']=row[str(year)]
    with (DATA/'observed_intervals.csv').open() as f:
        records=[]
        for r in csv.DictReader(f):
            if int(r['year'])!=year or r.get('candidate_interval','True')!='True':continue
            if r['site'] not in geom or 'aspect_deg' not in geom[r['site']]:raise ValueError('missing annual site geometry')
            r.update(geom[r['site']]);a,b=stamp(r['utc0']),stamp(r['utc1'])
            if a>=b or float(r['dx_horizontal_m'])<=0:continue
            records.append(dict(r,t0=a,t1=b,normal_m=float(r['dx_horizontal_m'])*math.sin(math.radians(float(r['slope_deg'])))))
        return records

def weather(year,frame=None):
    if frame is None:frame=pd.read_csv(RAW/'AWS'/f'Debris_EB_{year}.csv')
    times=pd.DatetimeIndex(pd.to_datetime(frame['Date_Time'],utc=True))
    required=['T_2m','RH_2m','Wspd','SWin','LWin']
    if not np.isfinite(frame[required].to_numpy(float)).all():raise ValueError('missing core weather')
    if ((frame.RH_2m<0)|(frame.RH_2m>100)|(frame.Wspd<0)|(frame.LWin<0)).any():raise ValueError('invalid physical weather values')
    ghi=np.maximum(frame.SWin.to_numpy(float),0.)
    solar=pvlib.solarposition.get_solarposition(times,61.4830,-142.9109,altitude=559.)
    split=pvlib.irradiance.erbs(ghi,solar.zenith,times)
    centers=times.asi8.astype(float)/1e9
    Ta=frame.T_2m.to_numpy(float);RH=frame.RH_2m.to_numpy(float)/100;wind=frame.Wspd.to_numpy(float)
    pressure=101325*(1-2.25577e-5*559.)**5.2559
    ea=RH*610.78*np.exp(17.27*Ta/(Ta+237.3));qs=0.622*610.78/(pressure-.378*610.78)
    qa=.622*ea/(pressure-.378*ea);rhoair=pressure/(287.05*(Ta+TM))
    turbulent=rhoair*wind*(1005.*Ta+2.5e6*(qa-qs))
    longwave=.97*(frame.LWin.to_numpy(float)-5.670374419e-8*TM**4)
    return dict(centers=centers,Ta=Ta,ghi=ghi,solar=solar,split=split,turbulent=turbulent,longwave=longwave,
                negative_SW_clipped_for_flux=int((frame.SWin<0).sum()))

def face(w,slope,aspect):
    x=pvlib.irradiance.get_total_irradiance(slope,aspect,w['solar'].zenith,w['solar'].azimuth,
                                         w['split']['dni'],w['ghi'],w['split']['dhi'],albedo=.2,model='isotropic')
    radiation=np.asarray(x['poa_global'],float)
    if not np.isfinite(radiation).all():raise ValueError('invalid radiation reconstruction')
    return radiation

def rate(w,sw,params):
    alpha,ch=params
    return np.maximum((1-alpha)*sw+w['longwave']+ch*w['turbulent'],0.)/(RHO*LATENT)

def main():
    protocol=json.loads((ROOT/'protocol.json').read_text())
    if protocol['training_year']!=2021 or protocol['test_year']!=2022:raise ValueError('unexpected protocol')
    train=intervals(2021);w1=weather(2021);features={};usable=[];excluded_train=[]
    dummy=Clock(w1['centers'],np.ones(len(w1['centers'])))
    for r in train:
        if r['t0']<dummy.edges[0] or r['t1']>dummy.edges[-1] or r['normal_m']<=0:
            excluded_train.append(dict(id=r['site']+':'+r['row0'],reason='no complete forcing or nonpositive retreat'));continue
        features.setdefault(r['site'],face(w1,float(r['slope_deg']),float(r['aspect_deg'])))
        usable.append(r)
    def residual(params):
        clocks={site:Clock(w1['centers'],rate(w1,sw,params)) for site,sw in features.items()}
        return np.array([clocks[r['site']].integral(r['t0'],r['t1'])-r['normal_m'] for r in usable])
    fits=[least_squares(residual,start,bounds=([.05,.0001],[.6,.02]),xtol=1e-11,ftol=1e-11,gtol=1e-11) for start in ([.2,.002],[.4,.005],[.1,.01])]
    fit=min(fits,key=lambda x:float(np.dot(x.fun,x.fun)))
    const=sum(r['normal_m'] for r in usable)/sum(r['t1']-r['t0'] for r in usable)
    degree=Clock(w1['centers'],np.maximum(w1['Ta'],0.)/86400)
    x=np.array([degree.integral(r['t0'],r['t1']) for r in usable]);y=np.array([r['normal_m'] for r in usable])
    pdd=float(np.dot(x,y)/np.dot(x,x))
    calibration=dict(protocol_sha256=hashlib.sha256((ROOT/'protocol.json').read_bytes()).hexdigest(),
                     training_year=2021,training_intervals=len(usable),excluded=excluded_train,
                     training_sites=sorted(features),albedo=float(fit.x[0]),neutral_exchange=float(fit.x[1]),
                     training_squared_error_m2=float(np.dot(fit.fun,fit.fun)),constant_normal_rate_m_s=const,
                     degree_factor_m_per_Kday=pdd,fit_success=bool(fit.success),
                     model_scope='atmosphere-only component baselines; unknown hydrological heat is not observed zero')
    # Save calibration BEFORE loading test-year target records or scoring.
    (ROOT/'CALIBRATION_2021.json').write_text(json.dumps(calibration,indent=2)+'\n')
    test=intervals(2022);w2=weather(2022);raw=[];exclusions=[];clocks={}
    support=Clock(w2['centers'],np.ones(len(w2['centers'])))
    for r in test:
        if r['site'] not in protocol['test_sites']:continue
        ident=r['site']+':'+r['row0']+'-'+r['row1']
        reason=None
        if abs(float(r['dx_horizontal_m'])-.1)>1e-9:reason='not an exact10cm marker interval'
        elif r['t0']<support.edges[0] or r['t1']>support.edges[-1]:reason='observed interval outside forcing support'
        if reason:exclusions.append(dict(id=ident,reason=reason));continue
        if r['site'] not in clocks:
            sw=face(w2,float(r['slope_deg']),float(r['aspect_deg']))
            clocks[r['site']]={'training_pooled_constant_normal_rate':Clock(w2['centers'],np.full(len(sw),const)),
                              'training_positive_degree_day':Clock(w2['centers'],pdd*np.maximum(w2['Ta'],0.)/86400),
                              'uncalibrated_neutral_energy_balance':Clock(w2['centers'],rate(w2,sw,(.2,.002))),
                              'training_calibrated_neutral_energy_balance':Clock(w2['centers'],rate(w2,sw,fit.x))}
        target=.1*math.sin(math.radians(float(r['slope_deg'])))
        observed=(r['t1']-r['t0'])/3600
        predicted={name:clock.arrival(r['t0'],target) for name,clock in clocks[r['site']].items()}
        raw.append(dict(id=ident,site=r['site'],source_rows=[r['row0'],r['row1']],start_utc=r['utc0'],observed_end_utc=r['utc1'],
                        observed_duration_h=observed,normal_target_m=target,
                        predicted_duration_h={k:(v-r['t0'])/3600 if v is not None else None for k,v in predicted.items()}))
    metrics={};models=protocol['models']
    for site in protocol['test_sites']+['pooled']:
        site_rows=[r for r in raw if site=='pooled' or r['site']==site]
        common=[r for r in site_rows if all(v is not None for v in r['predicted_duration_h'].values())]
        metrics[site]=dict(eligible=len(site_rows),paired_complete=len(common),models={})
        baseline=None
        for name in models:
            errors=np.array([r['predicted_duration_h'][name]-r['observed_duration_h'] for r in common])
            observed=np.array([r['observed_duration_h'] for r in common])
            if not len(errors):raise RuntimeError('no complete paired field events')
            mae=float(np.mean(abs(errors)))
            if baseline is None:baseline=mae
            metrics[site]['models'][name]=dict(MAE_h=mae,median_absolute_error_h=float(np.median(abs(errors))),
                RMSE_h=float(np.sqrt(np.mean(errors**2))),bias_h=float(errors.mean()),
                MAPE=float(np.mean(abs(errors)/observed)),fraction_within_2h=float(np.mean(abs(errors)<=2)),
                fraction_within20percent=float(np.mean(abs(errors)/observed<=.2)),
                MAE_skill_vs_training_constant=1-mae/baseline if baseline>0 else None,
                censored=sum(r['predicted_duration_h'][name] is None for r in site_rows))
    result=dict(evaluation='retrospective observed-weather hindcast, not an issued forecast',
                calibration=calibration,test_year=2022,metrics=metrics,exclusions=exclusions,events=raw,
                test_histories=2,claims='partial atmospheric component benchmark, not validated coupled model or professional-code reproduction',
                unknowns=['stream heat','terrain radiation/shading','changing slope/aspect within season','stake spatial continuity','hourly averaging convention'],
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                runtime_versions=dict(numpy=np.__version__,pandas=pd.__version__,pvlib=pvlib.__version__))
    (ROOT/'FIELD_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(calibration=calibration,metrics=metrics,exclusions=exclusions),indent=2))

if __name__=='__main__':main()
