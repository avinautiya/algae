"""Well-mixed reservoir + two ideal melting outlet controls.

Independent conditional reduced forecast, not spatial CFD or field validation.
Prescribed removal changes a uniform deposit; grain transport is not solved.
"""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

RHO_W=1000.
CP=4217.
RHO_I=917.
LATENT=334000.
G=9.81


def number(v, name, lower=None, positive=False):
    if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):
        raise ValueError(name+' must be finite numeric')
    if (lower is not None and v<lower) or (positive and v<=0):
        raise ValueError(name+' out of range')
    return float(v)


def validate(c):
    required={'pond_area_m2','pond_floor_m','initial_stage_m','initial_temperature_K',
              'horizon_s','output_step_s','routes','forcing'}
    allowed=required|{'removals','evidence','rtol','atol'}
    if not isinstance(c,dict) or set(c)-allowed or required-set(c):
        raise ValueError('unknown or missing top-level fields')
    A=number(c['pond_area_m2'],'pond_area',positive=True)
    floor=number(c['pond_floor_m'],'pond_floor')
    stage=number(c['initial_stage_m'],'initial_stage')
    if stage<=floor:raise ValueError('initial reservoir must contain water')
    number(c['initial_temperature_K'],'initial superheat',lower=0)
    horizon=number(c['horizon_s'],'horizon',positive=True)
    number(c['output_step_s'],'output_step',positive=True)
    if horizon/c['output_step_s']>1000000:raise ValueError('too many requested output samples')
    if not isinstance(c['routes'],list) or len(c['routes'])!=2:raise ValueError('exactly two routes required')
    keys={'ice_z_m','width_m','ice_area_m2','heat_transfer_W_m2_K','debris_thickness_m','debris_conductivity_W_m_K'}
    for route in c['routes']:
        if not isinstance(route,dict) or set(route)!=keys:raise ValueError('route keys invalid')
        if number(route['ice_z_m'],'ice_z')<=floor:raise ValueError('ice controls must be above fixed pond floor')
        for k in keys-{'ice_z_m'}:number(route[k],k,lower=0,positive=k in {'width_m','ice_area_m2','debris_conductivity_W_m_K'})
    if not isinstance(c['forcing'],list) or not c['forcing']:raise ValueError('forcing required')
    times=[]
    for f in c['forcing']:
        if not isinstance(f,dict) or set(f)!={'time_s','inflow_m3_s','temperature_K'}:raise ValueError('forcing keys invalid')
        times.append(number(f['time_s'],'forcing time',lower=0))
        number(f['inflow_m3_s'],'inflow',lower=0);number(f['temperature_K'],'inflow superheat',lower=0)
    if times[0]!=0 or any(a>=b for a,b in zip(times,times[1:])) or times[-1]>=horizon:
        raise ValueError('forcing starts at zero, strictly increases, and stays before horizon')
    removals=c.get('removals',[])
    if not isinstance(removals,list):raise ValueError('removals must be a list')
    for e in removals:
        if not isinstance(e,dict) or set(e)!={'time_s','route'}:raise ValueError('removal keys invalid')
        t=number(e['time_s'],'removal time',lower=0)
        if t>horizon or type(e['route']) is not int or e['route'] not in (0,1):raise ValueError('invalid removal')
    for key,default in [('rtol',1e-8),('atol',1e-10)]:number(c.get(key,default),key,positive=True)
    return c


def forecast(config):
    c=validate(config)
    A=c['pond_area_m2'];floor=c['pond_floor_m'];routes=c['routes']
    horizon=c['horizon_s'];step=c['output_step_s']
    volumes=(c['initial_stage_m']-floor)*A
    # V, reservoir sensible U, z1,z2, Vin,Vo1,Vo2,Ein,Eout,Elat1,Elat2.
    y=np.array([volumes,RHO_W*CP*volumes*c['initial_temperature_K'],
                routes[0]['ice_z_m'],routes[1]['ice_z_m'],0,0,0,0,0,0,0.],float)
    initial=y.copy();rows=[];work=0
    boundaries=sorted({0.,float(horizon)}|{float(f['time_s']) for f in c['forcing']}|
                      {float(e['time_s']) for e in c.get('removals',[])})

    def settings(t):
        forcing=max((f for f in c['forcing'] if f['time_s']<=t),key=lambda f:f['time_s'])
        thickness=[r['debris_thickness_m'] for r in routes]
        for event in c.get('removals',[]):
            if event['time_s']<=t:thickness[event['route']]=0.
        return forcing,thickness

    def rates(state,t):
        f,thickness=settings(t)
        V,U=state[:2];stage=floor+max(V,0)/A
        T=max(U,0)/(RHO_W*CP*V) if V>1e-15 else 0.
        flows=[];powers=[]
        for j,r in enumerate(routes):
            crest=state[2+j]+thickness[j]
            q=r['width_m']*math.sqrt(G)*(2/3)**1.5*max(stage-crest,0)**1.5
            h=r['heat_transfer_W_m2_K']
            conductance=0. if h==0 else 1/(1/h+thickness[j]/r['debris_conductivity_W_m_K'])
            p=RHO_W*CP*q*T*(-math.expm1(-conductance*r['ice_area_m2']/(RHO_W*CP*q))) if q>0 else 0.
            flows.append(q);powers.append(p)
        return f,stage,T,np.array(flows),np.array(powers)

    def emit(t,state):
        check_state(state)
        _,stage,T,q,p=rates(state,t)
        meltwater=p/(LATENT*RHO_W)
        downstream=q+meltwater
        downstream_T=np.divide(RHO_W*CP*q*T-p,RHO_W*CP*downstream,
                               out=np.zeros(2),where=downstream>0)
        rows.append(dict(time_s=float(t),stage_m=float(stage),reservoir_volume_m3=float(state[0]),
                         reservoir_sensible_energy_J=float(state[1]),
                         roundoff_clipping_applied=bool(state[0]<0 or state[1]<0),
                         reservoir_temperature_K=float(T),ice_z_m=state[2:4].tolist(),
                         route_flow_m3_s=q.tolist(),ice_power_W=p.tolist(),
                         downstream_route_flow_m3_s=downstream.tolist(),
                         downstream_route_temperature_K=downstream_T.tolist(),
                         toe_fraction=float(q[0]/q.sum()) if q.sum()>0 else None))

    def check_state(state):
        if not np.isfinite(state).all() or state[0]<-1e-10 or state[1]<-1e-6:
            raise RuntimeError('nonphysical reservoir state at accepted/output sample')

    for start,end in zip(boundaries,boundaries[1:]):
        f,thickness=settings(start)
        # Freeze boundary protocol within each exactly event-aligned segment.
        sample_t=start
        def rhs(t,state):
            _,_,T,q,p=rates(state,sample_t)
            dz=[-p[j]/(RHO_I*LATENT*routes[j]['ice_area_m2']) for j in range(2)]
            ein=RHO_W*CP*f['inflow_m3_s']*f['temperature_K']
            sensible_out=RHO_W*CP*q.sum()*T
            return np.array([f['inflow_m3_s']-q.sum(),ein-sensible_out,*dz,
                             f['inflow_m3_s'],q[0],q[1],ein,sensible_out-p.sum(),p[0],p[1]])
        def floor_event(t,state):return min(state[2:4])-floor
        floor_event.terminal=True;floor_event.direction=-1
        # Include a global regular grid and every protocol boundary.
        k0=math.ceil(start/step);k1=math.floor(end/step)
        times=sorted({float(start),float(end)}|{k*step for k in range(k0,k1+1) if start<=k*step<=end})
        result=solve_ivp(rhs,(start,end),y,method='Radau',rtol=c.get('rtol',1e-8),
                         atol=c.get('atol',1e-10),dense_output=True,events=floor_event)
        if not result.success:raise RuntimeError('integration failed: '+result.message)
        if len(result.t_events[0]):raise RuntimeError('control reached fixed pond floor; model scope exhausted')
        for state in result.y.T:check_state(state)
        work+=result.nfev
        for t in times:
            if rows and t==rows[-1]['time_s']:rows.pop()
            emit(t,result.sol(t))
        y=result.y[:,-1]
    check_state(y)
    water=y[0]-initial[0]-y[4]+y[5]+y[6]
    heat=y[1]-initial[1]+y[8]+y[9]+y[10]-y[7]
    latent=np.array([RHO_I*LATENT*r['ice_area_m2']*(initial[2+j]-y[2+j]) for j,r in enumerate(routes)])
    phase_errors=latent-y[9:11]
    heat_scale=max(1.,abs(initial[1])+abs(y[7]))
    water_scale=max(1.,abs(initial[0])+abs(y[4]))
    checks=dict(reservoir_water_closed=bool(abs(water)<=1e-8*water_scale),
                total_heat_closed=bool(abs(heat)<=1e-8*heat_scale),
                latent_exchange_closed=bool(np.max(np.abs(phase_errors))<=1e-8*heat_scale))
    if not all(checks.values()):raise RuntimeError('budget verification failed: '+str(checks))
    return dict(schema_version=1,input_sha256=hashlib.sha256(json.dumps(c,sort_keys=True,allow_nan=False).encode()).hexdigest(),
                model_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                evidence='conditional reduced-model forecast; not spatial CFD or externally validated',
                scope=['well-mixed fixed-hypsometry reservoir','ideal broad-crested parallel controls',
                       'quasi-steady plug-flow heat loss','uniform prescribed debris cover/removal',
                       'meltwater joins downstream and does not feed back into reservoir'],
                series=rows,rhs_evaluations=work,checks=checks,
                numerics=dict(method='Radau',rtol=c.get('rtol',1e-8),scalar_atol=c.get('atol',1e-10),
                              note='mixed SI units; refine solution observables, not just algebraic ledgers'),
                applicability=dict(meltwater_to_route_inflow_ratio_upper_bound=CP*max(
                    c['initial_temperature_K'],max(f['temperature_K'] for f in c['forcing']))/LATENT,
                    ratio_scope='passive well-mixed maximum principle and all advected heat assigned to latent melt; not a hydraulic-error estimate'),
                ledger=dict(reservoir_water_residual_m3=float(water),total_heat_residual_J=float(heat),
                            route_latent_residual_J=phase_errors.tolist(),
                            inlet_water_m3=float(y[4]),outlet_water_from_reservoir_m3=y[5:7].tolist(),
                            downstream_meltwater_m3=(y[9:11]/LATENT/RHO_W).tolist(),
                            ice_latent_energy_J=y[9:11].tolist(),inlet_sensible_energy_J=float(y[7]),
                            downstream_sensible_energy_J=float(y[8]),final_reservoir_sensible_energy_J=float(y[1])))
