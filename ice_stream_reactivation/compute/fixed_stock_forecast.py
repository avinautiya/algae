"""One-quadrature forecast for fixed, unexported debris after fines clear.

Same reduced closures; not a new fundamental instability theorem. Invalid if
coarse grains move, forcing changes, or spatial water paths defeat mean crest.
"""
import math
import numpy as np
from scipy.integrate import quad
from model import partition


def predict_loss(c, solid_stock_m3, toe_ice_z=None, bypass_ice_z=None):
    if c.allow_supply:raise ValueError('fixed-stock forecast requires no incoming debris')
    S=np.asarray(solid_stock_m3,dtype=float);D=np.asarray(c.grain_sizes)
    if S.shape!=D.shape or np.any(S<0) or not np.all(np.isfinite(S)):
        raise ValueError('nonnegative finite grain stocks required')
    A=c.width*c.length;As=(1-c.porosity)*A
    zt=c.initial_toe_ice_z if toe_ice_z is None else toe_ice_z
    zb=c.initial_bypass_z if bypass_ice_z is None else bypass_ice_z
    K=(2/3)**1.5*math.sqrt(c.gravity)*c.width
    target=(c.discharge/K)**(2/3)
    if not math.isfinite(zt+zb):raise ValueError('finite ice elevations required')
    g0=zt-zb+S.sum()/As
    footprint=S/D;total=footprint.sum();cover=min(total/A,1.)
    if c.h_water:
        covered=float(np.dot(footprint/total,1/(1/c.h_water+D/c.debris_k))) if total else 0.
        ht=(1-cover)*c.h_water+cover*covered
    else:ht=0.
    if not c.footprint_shield:
        ht=1/(1/c.h_water+S.sum()/As/c.debris_k) if c.h_water else 0.
    if not c.allow_shield:ht=c.h_water
    def power(q,h):
        cpq=c.water_density*c.water_cp*q
        return cpq*c.water_temperature*(-math.expm1(-h*A/cpq)) if cpq>0 else 0.
    def drift(gap):
        qt,qb,_=partition(c.discharge,gap,0.,K,K)
        return (power(qb,c.h_water)-power(qt,ht)-c.extra_toe_heat_W_m2*A)/(c.ice_density*c.latent*A)
    if g0>=target:return dict(status='initially dry',loss_time_days=0.)
    if not c.allow_ice_incision:return dict(status='no ice-bed evolution',loss_time_days=None)
    r0=drift(g0)
    if r0<=0:return dict(status='no forward dry-threshold crossing in the fixed-stock autonomous model',loss_time_days=None)
    r1=drift(target)
    seconds,err=quad(lambda gap:1/drift(gap),g0,target,epsabs=.1,epsrel=1e-8,limit=120)
    return dict(status='finite loss in fixed-stock autonomous model',loss_time_days=seconds/86400,
                mathematical_time_bracket_days=[(target-g0)/r1/86400,(target-g0)/r0/86400],
                quadrature_error_estimate_days=err/86400,initial_coverage=float(cover),
                scope='fixed stock, heat transfer, width, forcing and mean-crest hydraulic model')

if __name__=='__main__':
    from pathlib import Path
    import json
    from model import Case
    rows=[]
    for frac in (.2,.4,.6,.8,1.):
        c=Case(discharge=.07,initial_debris_depth=.2,allow_supply=False,source_fractions=((1-frac)/2,(1-frac)/2,frac))
        total=c.initial_debris_depth*(1-c.porosity)*c.width*c.length
        result=predict_loss(c,[0.,0.,total*frac])
        rows.append(dict(coarse_fraction=frac,forecast=result));print(rows[-1],flush=True)
    (Path(__file__).parent/'records/fixed_stock_forecast.json').write_text(json.dumps(rows,indent=2))
