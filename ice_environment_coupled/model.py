"""Coupled reduced cliff forecast with conservative environmental energy ports.

Not spatial CFD: a thin upstream obstruction feeds a bare downstream floor;
the upper face may carry DRY debris. No under-debris floor thermal shortcut.
"""
from dataclasses import asdict,replace
import hashlib,json,math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).parent/'vendor'))
from stepped_cliff import SteppedCliff,advance_latent
from transmission import Route,solve_common_stage,hydrostatic_diagnostics
from heat_delivery import water_heat_rates
from energy_ports.energy import State,Parameters,AreaContent,remap_area,step,step_bare

RHO_I=917.;RHO_W=1000.;L=334000.;CP=4217.

def contact_depth(q,width,n,friction_slope):
    """Rectangular Manning rating for an EXPLICIT imposed friction slope.

    This is a downstream contact-depth closure, not the upstream pool stage,
    a resolved water profile or a slope inferred from the ice geometry.
    """
    if q<=0:return 0.
    def rating(h):return width*h*(width*h/(width+2*h))**(2/3)*math.sqrt(friction_slope)/n
    lo=0.;hi=1.
    while rating(hi)<q:hi*=2
    for _ in range(80):
        mid=(lo+hi)/2
        if rating(mid)<q:lo=mid
        else:hi=mid
    return (lo+hi)/2

def finite(v,name,low=None):
    if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or (low is not None and v<low):
        raise ValueError('invalid '+name)
    return float(v)

def forecast(c):
    required={'geometry','physics','periods','step_s','removal_time_s','provenance'}
    if not isinstance(c,dict) or set(c)!=required:raise ValueError('unknown/missing configuration fields')
    if not isinstance(c['provenance'],str) or not c['provenance'].strip():raise ValueError('input provenance required')
    g=SteppedCliff(**c['geometry']).validate();p=c['physics']
    keys={'bypass_floor_m','bypass_floor_limit_m','bypass_width_m','bypass_area_m2','weir_coefficient_sqrt_m_s',
          'obstruction_height_m','aperture_m','discharge_coefficient','darcy_conductance_m2_s',
          'water_h_W_m2_K','face_cover_fraction','debris_capacity_J_m2_K','debris_conductance_W_m2_K',
          'initial_debris_temperature_C','new_debris_temperature_C','initial_surface_cold_J_m2',
          'contact_manning_n','contact_friction_slope',
          'new_surface_cold_J_m2','ice_subcooling_K','ice_heat_capacity_J_kg_K','max_meltwater_fraction'}
    if not isinstance(p,dict) or set(p)!=keys:raise ValueError('unknown/missing physics fields')
    signed={'bypass_floor_m','bypass_floor_limit_m','initial_debris_temperature_C','new_debris_temperature_C'}
    for k,v in p.items():finite(v,k,None if k in signed else 0)
    for k in ('bypass_width_m','bypass_area_m2','weir_coefficient_sqrt_m_s','debris_capacity_J_m2_K','ice_heat_capacity_J_kg_K','contact_manning_n','contact_friction_slope'):
        if p[k]<=0:raise ValueError(k+' must be positive')
    f=p['face_cover_fraction']
    if not 0<=f<1:raise ValueError('face cover fraction must be >=0 and <1')
    if p['bypass_floor_limit_m']>=p['bypass_floor_m']:raise ValueError('bypass domain invalid')
    dtmax=finite(c['step_s'],'step_s',0)
    if dtmax==0:raise ValueError('positive step required')
    if not isinstance(c['periods'],list) or not c['periods']:raise ValueError('forcing periods required')
    end=0.;periods=[]
    for row in c['periods']:
        if set(row)!={'duration_s','inflow_m3_s','water_superheat_K','face_atmospheric_net_W_m2',
                     'floor_atmospheric_net_W_m2','bypass_atmospheric_net_W_m2'}:raise ValueError('forcing schema')
        for k,v in row.items():finite(v,k,0 if k in ('duration_s','inflow_m3_s','water_superheat_K') else None)
        if row['duration_s']==0:raise ValueError('positive forcing duration required')
        periods.append((end,end+row['duration_s'],row));end+=row['duration_s']
    removal=c['removal_time_s']
    if removal is not None:
        finite(removal,'removal time',0)
        if removal>end:raise ValueError('removal beyond horizon')
    if math.ceil(end/dtmax)>200000:raise ValueError('too many steps')
    par=Parameters(p['debris_capacity_J_m2_K'],p['debris_conductance_W_m2_K'],L)
    areas={'face_bare':g.face_area_m2*(1-f),'floor':g.floor_area_m2,'bypass':p['bypass_area_m2']}
    if f>0:areas['face_debris']=g.face_area_m2*f
    states={k:State(a,par.capacity_j_m2_k*a*p['initial_debris_temperature_C'] if k=='face_debris' else 0.,
                    a*p['initial_surface_cold_J_m2']) for k,a in areas.items()}
    def stored():return sum(s.debris_energy_j-s.cold_debt_j for s in states.values())
    initial_stored=stored();initial_ice=g.ice_volume_m3;zb=p['bypass_floor_m'];zb0=zb
    totals=dict(atmospheric_J=0.,water_in_J=0.,water_out_J=0.,latent_J=0.,fresh_ice_warming_J=0.,
                remap_enthalpy_J=0.,face_melt_m3=0.,floor_melt_m3=0.,bypass_melt_m3=0.,water_in_m3=0.,
                water_out_m3=0.,phase_water_m3=0.)
    cold_vol=RHO_I*p['ice_heat_capacity_J_kg_K']*p['ice_subcooling_K'];rows=[];t=0.
    def phase(k,res):
        volume=res.latent_j/(RHO_I*L+cold_vol)
        latent=volume*RHO_I*L
        states[k]=replace(res.state,cumulative_latent_j=states[k].cumulative_latent_j+latent)
        totals['latent_J']+=latent;totals['fresh_ice_warming_J']+=volume*cold_vol
        return latent,volume
    for start,stop,forcing in periods:
        while t<stop:
            dt=min(dtmax,stop-t)
            if removal is not None and t<removal<t+dt:dt=removal-t
            if t+dt<=t:raise ValueError('forcing interval below floating-point time resolution')
            removed=removal is not None and t>=removal
            D=0. if removed else p['obstruction_height_m']
            aperture=0. if removed else p['aperture_m']
            routes=[Route('toe',g.floor_z_m,D,g.extrusion_width_m,p['weir_coefficient_sqrt_m_s'],aperture,
                          p['discharge_coefficient'],0. if removed else p['darcy_conductance_m2_s']),
                    Route('bypass',zb,0.,p['bypass_width_m'],p['weir_coefficient_sqrt_m_s'])]
            qin=forcing['inflow_m3_s'];Tw=forcing['water_superheat_K']
            hydro=solve_common_stage(routes,qin,min(g.floor_z_m,zb)-1.,flow_tolerance_m3_s=1e-12)
            qt=hydro.route_flows_m3_s['toe'];qb=hydro.route_flows_m3_s['bypass']
            depth=contact_depth(qt,g.extrusion_width_m,p['contact_manning_n'],p['contact_friction_slope'])
            contact_stage=g.floor_z_m+depth
            velocity=qt/(g.extrusion_width_m*depth) if qt>1e-14 and depth>0 else 0.
            required_head=contact_stage+velocity**2/(2*9.81)+p['contact_friction_slope']*(g.face_x_m-g.upstream_edge_x_m)
            if qt>1e-14 and required_head>hydro.stage_m+1e-10:
                raise ValueError('route/contact closures require more head than the upstream pool supplies')
            indicators=hydrostatic_diagnostics(depth,g.face_x_m-g.upstream_edge_x_m,
                       qt/(g.extrusion_width_m*depth),0.,None) if qt>1e-14 and depth>0 else None
            cover_bottom=g.crest_z_m-f*(g.crest_z_m-g.floor_z_m)
            wet=max(0.,min(cover_bottom,contact_stage)-g.floor_z_m)*g.extrusion_width_m if qt>1e-14 else 0.
            if qt>1e-14 and f>0 and contact_stage>cover_bottom:
                raise ValueError('submerged face debris requires a separate thermal boundary closure')
            heat=water_heat_rates(qt,Tw,p['water_h_W_m2_K']*g.floor_area_m2,
                                  p['water_h_W_m2_K']*wet)
            byheat=water_heat_rates(qb,Tw,p['water_h_W_m2_K']*p['bypass_area_m2'],0.)
            at={k:0. for k in states};at['face_bare']=forcing['face_atmospheric_net_W_m2']*(states['face_bare'].area_m2-wet)
            at['floor']=forcing['floor_atmospheric_net_W_m2']*g.floor_area_m2 if qt<=1e-14 else 0.
            at['bypass']=forcing['bypass_atmospheric_net_W_m2']*p['bypass_area_m2'] if qb<=1e-14 else 0.
            if f>0:at['face_debris']=forcing['face_atmospheric_net_W_m2']*states['face_debris'].area_m2
            water={'floor':heat['floor_water_heat_W'],'face_bare':heat['face_water_heat_W'],
                   'bypass':byheat['floor_water_heat_W'],'face_debris':0.}
            latent={};vol={}
            for k,state in list(states.items()):
                res=step(state,par,dt,at[k]/state.area_m2,0.) if k=='face_debris' else step_bare(state,dt,at[k]/state.area_m2,water[k]/state.area_m2)
                if abs(res.residual_j)>1e-7*max(1.,abs(res.atmospheric_j)+abs(res.water_j)):raise RuntimeError('port energy failure')
                latent[k],vol[k]=phase(k,res)
            face_latent=latent['face_bare']+latent.get('face_debris',0.)
            g,geom=advance_latent(g,face_latent,latent['floor'])
            zb-=vol['bypass']/p['bypass_area_m2']
            if zb<=p['bypass_floor_limit_m']:raise ValueError('bypass ice domain exhausted')
            phasewater=RHO_I*sum(vol.values())/RHO_W
            ratio=phasewater/(qin*dt) if qin>0 else (float('inf') if phasewater>0 else 0.)
            if ratio>p['max_meltwater_fraction']:raise ValueError('meltwater hydraulic feedback exceeds declared reduced-model scope')
            newareas={'face_bare':g.face_area_m2*(1-f),'floor':g.floor_area_m2,'bypass':p['bypass_area_m2']}
            if f>0:newareas['face_debris']=g.face_area_m2*f
            for k,s in list(states.items()):
                if newareas[k]==s.area_m2:continue
                if newareas[k]<s.area_m2:raise RuntimeError('unimplemented material contraction')
                content=AreaContent(par.capacity_j_m2_k*p['new_debris_temperature_C'] if k=='face_debris' else 0.,
                                    p['new_surface_cold_J_m2'],0.)
                rem=remap_area(s,newareas[k],incoming=content);states[k]=rem.state
                totals['remap_enthalpy_J']+=rem.debris_transfer_j-rem.cold_debt_transfer_j
            totals['atmospheric_J']+=sum(at.values())*dt
            totals['water_in_J']+=RHO_W*CP*(qt+qb)*Tw*dt
            totals['water_out_J']+=(heat['outgoing_water_heat_W']+byheat['outgoing_water_heat_W'])*dt
            totals['water_in_m3']+=qin*dt
            totals['water_out_m3']+=(qt+qb)*dt+phasewater
            totals['phase_water_m3']+=phasewater
            totals['face_melt_m3']+=vol['face_bare']+vol.get('face_debris',0.)
            totals['floor_melt_m3']+=vol['floor'];totals['bypass_melt_m3']+=vol['bypass']
            residual=(totals['atmospheric_J']+totals['water_in_J']-totals['water_out_J']+
                      totals['remap_enthalpy_J']-(stored()-initial_stored)-totals['latent_J']-totals['fresh_ice_warming_J'])
            t+=dt
            rows.append(dict(time_s=t,stage_m=hydro.stage_m,toe_flow_m3_s=qt,bypass_flow_m3_s=qb,
                             upstream_stage_m=hydro.stage_m,contact_stage_m=contact_stage,contact_depth_m=depth,
                             contact_required_head_m=required_head,
                             hydraulic_sample_time_s=t-dt,hydrostatic_indicators=indicators,
                             toe_fraction=qt/qin if qin>0 else None,removed=removed,geometry=asdict(g),bypass_floor_m=zb,
                             wet_bare_face_area_m2=wet,face_melt_m3=totals['face_melt_m3'],
                             floor_melt_m3=totals['floor_melt_m3'],bypass_melt_m3=totals['bypass_melt_m3'],
                             global_energy_residual_J=residual,routing_residual_m3_s=hydro.residual_m3_s,
                             meltwater_to_inflow_ratio=ratio,port_states={k:asdict(v) for k,v in states.items()}))
    phase_geom=initial_ice-g.ice_volume_m3+(zb0-zb)*p['bypass_area_m2']
    volumes=totals['face_melt_m3']+totals['floor_melt_m3']+totals['bypass_melt_m3']
    scale=max(1.,abs(totals['atmospheric_J'])+totals['water_in_J'])
    if abs(residual)>1e-9*scale or abs(phase_geom-volumes)>1e-9*max(1.,volumes):raise RuntimeError('global conservation failure')
    return dict(schema_version=1,input_sha256=hashlib.sha256(json.dumps(c,sort_keys=True,allow_nan=False).encode()).hexdigest(),
                evidence='conditional coupled reduced-model forecast; not spatial CFD or field validation',
                provenance=c['provenance'],series=rows,totals=totals,
                checks=dict(global_energy_residual_J=residual,geometry_phase_volume_residual_m3=phase_geom-volumes,
                            water_bookkeeping_residual_m3=totals['water_out_m3']-totals['water_in_m3']-totals['phase_water_m3']),
                scope=['thin effective upstream obstacle; bare downstream floor',
                       'upper dry face debris; submerged debris rejected',
                       'quasi-steady common stage and uniform planar retreat',
                       'downstream face wetting from imposed-friction-slope rectangular Manning rating, not upstream stage',
                       'prescribed lumped face fractions; no tracked bare/debris boundary advection',
                       'explicit incoming surface reservoirs and bulk fresh-ice sensible heat',
                       'meltwater leaves downstream; small-feedback bound enforced',
                       'prescribed net atmospheric face flux; no inferred runoff or water temperature',
                       'no collapse, glacier advection or resolved cold-ice PDE'])

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.input.resolve()==args.output.resolve():parser.error('output cannot replace input')
    result=forecast(json.loads(args.input.read_text()));args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result['checks'],indent=2))
