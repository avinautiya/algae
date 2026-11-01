"""Independent adaptive time integration + Brent hydraulics of finite-stock case.

Shares declared constitutive parameters, not production stepping or partition.
Numerical verification of reduced equations, not independent physical validation.
"""
import math,json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from model import Case


def reference(c,days=30,pulse=None):
    if c.extra_toe_heat_W_m2 or len(c.grain_sizes)!=3:raise ValueError('reference only implements zero extra heat and three size classes')
    if c.allow_supply:raise ValueError('reference only implements finite initial stock')
    A=c.width*c.length;As=(1-c.porosity)*A
    K=(2/3)**1.5*math.sqrt(c.gravity)*c.width
    D=np.array(c.grain_sizes);latent=c.ice_density*c.latent
    y=np.r_[c.initial_toe_ice_z,c.initial_bypass_z,c.initial_debris_depth*As*np.array(c.source_fractions),0.,0.]
    points=[0.,days*86400]
    if pulse:points+= [pulse['start_day']*86400,(pulse['start_day']+pulse['duration_day'])*86400]
    points=sorted(t for t in set(points) if 0<=t<=days*86400)
    events=[];nfev=0;minstock=0.
    for ta,tb in zip(points[:-1],points[1:]):
        Q=c.discharge
        if pulse and pulse['start_day']*86400 <=(ta+tb)/2<(pulse['start_day']+pulse['duration_day'])*86400:Q*=pulse['multiplier']
        def rhs(t,y):
            S=y[2:5];d=max(float(S.sum()/As),0.)
            zt=y[0]+d;zb=y[1]
            def balance(eta):return K*(max(eta-zt,0)**1.5+max(eta-zb,0)**1.5)-Q
            eta=brentq(balance,min(zt,zb),max(zt,zb)+(Q/K)**(2/3),xtol=1e-13)
            qt=K*max(eta-zt,0)**1.5;qb=K*max(eta-zb,0)**1.5
            depth=2*max(eta-zt,0)/3
            u=qt/(c.width*depth) if depth>0 else 0.
            def P(q,h):
                cpq=c.water_density*c.water_cp*q
                return cpq*c.water_temperature*(-math.expm1(-h*A/cpq)) if cpq>0 and h>0 else 0.
            ht=1/(1/c.h_water+d/c.debris_k) if c.h_water else 0.
            if c.footprint_shield and c.h_water:
                footprints=np.maximum(S,0)/D
                footprint=float(footprints.sum());cover=min(footprint/A,1.)
                covered=float(np.dot(footprints/footprint,1/(1/c.h_water+D/c.debris_k))) if footprint else 0.
                ht=(1-cover)*c.h_water+cover*covered
            Pt=P(qt,ht);Pb=P(qb,c.h_water)
            if not c.allow_ice_incision:Pt=Pb=0.
            mt,mb=Pt/(latent*A),Pb/(latent*A)
            sub=np.minimum(depth,D)
            normal=c.gravity*D**2*(c.grain_density*D-c.water_density*sub)
            U2=np.divide(2*c.friction*normal,c.water_density*c.drag*D*sub,out=np.full_like(D,np.inf),where=sub>0)
            speed=np.maximum(u-np.sqrt(np.maximum(U2,0)),0.)
            active=np.minimum(1.,D/max(d,1e-300))
            export=c.export_efficiency*active*speed/c.length*S if c.allow_export else np.zeros(3)
            return np.r_[-mt,-mb,-export,mt*A,mb*A]
        def dry(t,y):return y[0]+y[2:5].sum()/As-y[1]-(Q/K)**(2/3)
        sol=solve_ivp(rhs,(ta,tb),y,method='Radau',rtol=2e-9,
                      atol=np.array([1e-10,1e-10,1e-12,1e-12,1e-12,1e-10,1e-10]),
                      max_step=3600.,events=dry)
        if not sol.success:raise RuntimeError(sol.message)
        minstock=min(minstock,float(np.min(sol.y[2:5])))
        nfev+=sol.nfev;events.extend((sol.t_events[0]/86400).tolist());y=sol.y[:,-1]
    return dict(final_state=y.tolist(),toe_melt_m3=float(y[5]),bypass_melt_m3=float(y[6]),
                events_days=events,minimum_numerical_stock_m3=minstock,nfev=nfev)

if __name__=='__main__':
    rows=[];base=Case(footprint_shield=False,discharge=.07,initial_debris_depth=.2,allow_supply=False)
    for pulse in (None,dict(start_day=1.,multiplier=2.,duration_day=.25),dict(start_day=14.,multiplier=2.,duration_day=.25)):
        r=reference(base,pulse=pulse);rows.append(dict(pulse=pulse,result=r));print(rows[-1],flush=True)
    (Path(__file__).parent/'records/independent_reference.json').write_text(json.dumps(rows,indent=2))
