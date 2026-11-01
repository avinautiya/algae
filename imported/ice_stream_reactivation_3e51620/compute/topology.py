"""Ideal control-network experiment at equal debris area and solid volume.

Parallel longitudinal patches versus a transverse controlling covered reach.
Fixed immobile patch before prescribed complete removal: no transport claim.
Not resolved 2-D CFD. Transverse downstream bare reach is assumed noncontrolling.
"""
import math
import numpy as np
from model import Case


def run(c, layout, cover_fraction=.48, grain_depth=.2, removal_day=None, days=45., dt=225.):
    if layout not in ('corridor','spanning_strip') or not 0<cover_fraction<1:
        raise ValueError('invalid layout')
    if c.extra_toe_heat_W_m2 or min(grain_depth,days,dt)<=0:raise ValueError('unsupported forcing or invalid scale')
    if layout=='corridor':
        widths=np.array([(1-cover_fraction)*c.width,cover_fraction*c.width,c.width])
        lengths=np.full(3,c.length)
        bed=np.array([c.initial_toe_ice_z,c.initial_toe_ice_z,c.initial_bypass_z],dtype=float)
        covered=np.array([False,True,False]);toe=np.array([True,True,False])
    else:
        widths=np.array([c.width,c.width]);lengths=np.array([cover_fraction*c.length,c.length])
        bed=np.array([c.initial_toe_ice_z,c.initial_bypass_z],dtype=float)
        covered=np.array([True,False]);toe=np.array([True,False])
    areas=widths*lengths;K=(2/3)**1.5*math.sqrt(c.gravity)*widths
    melt=np.zeros_like(bed);wettime=0.;history=[];Q=c.discharge;latent=c.ice_density*c.latent
    heatin=0.;heatout=0.;minimum_bare_flow=float('inf')
    for step in range(math.ceil(days*86400/dt)):
        t=step*dt;h=min(dt,days*86400-t)
        present=removal_day is None or t<removal_day*86400
        obstacle=covered.astype(float)*grain_depth*present
        crest=bed+obstacle
        lo=float(np.min(crest));hi=float(np.max(crest)+(Q/np.min(K))**(2/3))
        for _ in range(44):
            eta=(lo+hi)/2
            if float(np.sum(K*np.maximum(eta-crest,0)**1.5))<Q:lo=eta
            else:hi=eta
        eta=(lo+hi)/2;q=K*np.maximum(eta-crest,0)**1.5
        ht=np.full_like(bed,c.h_water)
        if c.h_water:ht=np.where(covered & present,1/(1/c.h_water+grain_depth/c.debris_k),ht)
        capacity=c.water_density*c.water_cp*q
        ntu=np.divide(ht*areas,capacity,out=np.zeros_like(q),where=capacity>0)
        P=capacity*c.water_temperature*(-np.expm1(-ntu))
        if not c.allow_ice_incision:P*=0.
        m=P/(latent*areas)
        wettime+=h*(q[toe].sum()>1e-9*Q)
        if layout=='corridor':minimum_bare_flow=min(minimum_bare_flow,float(q[0]))
        if step%max(1,int(21600/dt))==0:
            history.append(dict(day=t/86400,flow=q.tolist(),ice_bed=bed.tolist(),crest=crest.tolist()))
        bed-=h*m;melt+=h*m*areas
        heatin+=h*c.water_density*c.water_cp*Q*c.water_temperature
        heatout+=h*(c.water_density*c.water_cp*Q*c.water_temperature-P.sum())
    energy_residual=heatin-heatout-latent*melt.sum()
    return dict(layout=layout,cover_fraction=cover_fraction,grain_depth=grain_depth,removal_day=removal_day,
                covered_area_m2=cover_fraction*c.width*c.length,
                solid_patch_volume_m3=grain_depth*cover_fraction*c.width*c.length,
                wet_days=wettime/86400,toe_melt_m3=float(melt[toe].sum()),
                final_ice_bed=bed.tolist(),minimum_bare_flow_m3_s=minimum_bare_flow if layout=='corridor' else None,
                energy_relative_residual=float(energy_residual/max(heatin,1.)),history=history,
                scope='ideal control-network topology comparison; downstream bare reach excluded from strip control melt')

if __name__=='__main__':
    from pathlib import Path
    import json
    c=Case(discharge=.07,allow_supply=False)
    rows=[]
    for f in (.2,.48,.8):
        for layout in ('corridor','spanning_strip'):
            for removal in (None,1.,14.):
                r=run(c,layout,cover_fraction=f,removal_day=removal)
                rows.append(r)
                print(f,layout,removal,r['wet_days'],r['toe_melt_m3'],flush=True)
    (Path(__file__).parent/'records/topology_comparison.json').write_text(json.dumps(rows,indent=2))
