"""Reduced cliff forecast with conservative thermal and geometric state.

Optional normal-resolved wet/dry debris and height-band face retreat. A thin
upstream obstruction feeds a bare downstream floor. Prescribed hydrostatic
closures and adiabatic ledges are not resolved free-surface flow or fracture.
"""
from dataclasses import asdict,replace
import hashlib,json,math
from pathlib import Path
import sys
from types import SimpleNamespace
if __package__:
    from .vendor.stepped_cliff import SteppedCliff,advance_latent
    from .vendor.transmission import Route,solve_common_stage,hydrostatic_diagnostics
    from .vendor.heat_delivery import water_heat_rates
    from .energy_ports.energy import State,Parameters,AreaContent,remap_area,step,step_bare
    from .thermal.wet_debris import advance as wet_advance, ThermalState, ThermalParameters
    from .transfer.local_transfer import evaluate as local_transfer, config_from_dict
    from .thermal.column_partition import ColumnPartition,remap as remap_columns
    from .face_thermal import advance_columns
    from .face_geometry.integration_adapter import initialize_profile,advance_with_partition
    from .face_geometry.face_profile import Profile,Band,face_partition as class_partition
    from .face_geometry.mesh import split_at_elevations
else:
    sys.path.insert(0,str(Path(__file__).parent/'vendor'))
    from stepped_cliff import SteppedCliff,advance_latent
    from transmission import Route,solve_common_stage,hydrostatic_diagnostics
    from heat_delivery import water_heat_rates
    from energy_ports.energy import State,Parameters,AreaContent,remap_area,step,step_bare
    from thermal.wet_debris import advance as wet_advance, ThermalState, ThermalParameters
    from transfer.local_transfer import evaluate as local_transfer, config_from_dict
    from thermal.column_partition import ColumnPartition,remap as remap_columns
    from face_thermal import advance_columns
    from face_geometry.integration_adapter import initialize_profile,advance_with_partition
    from face_geometry.face_profile import Profile,Band,face_partition as class_partition
    from face_geometry.mesh import split_at_elevations

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
    if not isinstance(c,dict) or set(c)-{'upgrades'}!=required:raise ValueError('unknown/missing configuration fields')
    upgrades=c.get('upgrades',{})
    if not isinstance(upgrades,dict) or set(upgrades)-{'wet_debris','face_cover_location','toe_heat_transfer','bypass_heat_transfer','bypass_friction_slope','meltwater_policy','debris_cells','debris_substep_fraction','covered_heat_policy','face_profile_bands','ledge_heat_policy','water_debris_coupling','initial_face_bands','profile_provenance','fixed_cover_interval_m','cover_interval_provenance','profile_alignment','max_profile_bands','profile_snapshot_stride','profile_ledger_detail'}:
        raise ValueError('unknown upgrade fields')
    cover_location=upgrades.get('face_cover_location','upper')
    if cover_location not in ('upper','lower'):raise ValueError('face cover location must be upper or lower')
    wet_mode=upgrades.get('wet_debris','reject')
    if wet_mode not in ('reject','lumped_temperate','resolved_temperate'):raise ValueError('unsupported wet debris mode')
    water_coupling=upgrades.get('water_debris_coupling','coupled')
    if water_coupling not in ('coupled','disabled_control'):raise ValueError('unsupported water/debris coupling')
    if water_coupling=='disabled_control' and wet_mode=='reject':raise ValueError('thermal intervention requires an explicit debris model')
    heat_policy=upgrades.get('covered_heat_policy','dry_only')
    if heat_policy not in ('dry_only','prescribed_all_cover'):raise ValueError('unsupported covered heat policy')
    if heat_policy!='dry_only' and wet_mode=='reject':raise ValueError('all-cover heat requires an explicit debris thermal model')
    debris_cells=upgrades.get('debris_cells',16)
    if isinstance(debris_cells,bool) or not isinstance(debris_cells,int) or not 4<=debris_cells<=256:raise ValueError('debris cells must be an integer from 4 to 256')
    substep_fraction=finite(upgrades.get('debris_substep_fraction',.01),'debris substep fraction',0)
    if not 0<substep_fraction<=.1:raise ValueError('debris substep fraction must be in (0,0.1]')
    profile_bands=upgrades.get('face_profile_bands')
    profile_alignment=upgrades.get('profile_alignment','surface_classes')
    if profile_alignment not in ('fixed','surface_classes'):raise ValueError('unknown profile alignment policy')
    snapshot_stride=upgrades.get('profile_snapshot_stride',1)
    if isinstance(snapshot_stride,bool) or not isinstance(snapshot_stride,int) or snapshot_stride<1:
        raise ValueError('profile snapshot stride must be a positive integer')
    ledger_detail=upgrades.get('profile_ledger_detail','full')
    if ledger_detail not in ('full','summary'):raise ValueError('unknown profile ledger detail')
    max_profile_bands=upgrades.get('max_profile_bands',4096)
    if isinstance(max_profile_bands,bool) or not isinstance(max_profile_bands,int) or not 4<=max_profile_bands<=16384:
        raise ValueError('maximum profile bands must be an integer in [4,16384]')
    if profile_bands is not None:
        if isinstance(profile_bands,bool) or not isinstance(profile_bands,int) or not 1<=profile_bands<=512:
            raise ValueError('profile bands must be an integer in [1,512]')
        if upgrades.get('ledge_heat_policy')!='adiabatic':raise ValueError('height profile requires explicit adiabatic ledge heat policy')
    elif 'ledge_heat_policy' in upgrades:raise ValueError('ledge policy supplied without a height profile')
    melt_policy=upgrades.get('meltwater_policy')
    if melt_policy is not None:
        if not isinstance(melt_policy,dict) or set(melt_policy)!={'mode','max_generated_discharge_m3_s','provenance'} or melt_policy['mode']!='downstream_free_drainage':
            raise ValueError('unsupported meltwater policy')
        if finite(melt_policy['max_generated_discharge_m3_s'],'meltwater drainage capacity',0)<=0:raise ValueError('positive drainage capacity required')
        if not isinstance(melt_policy['provenance'],str) or not melt_policy['provenance'].strip():raise ValueError('drainage assumption provenance required')
    if not isinstance(c['provenance'],str) or not c['provenance'].strip():raise ValueError('input provenance required')
    g=SteppedCliff(**c['geometry']).validate();p=c['physics']
    face_profile=initialize_profile(asdict(g),profile_bands) if profile_bands is not None else None
    initial_bands=upgrades.get('initial_face_bands')
    if initial_bands is not None:
        if face_profile is None:raise ValueError('initial bands require an explicit height-profile mode')
        if not isinstance(upgrades.get('profile_provenance'),str) or not upgrades['profile_provenance'].strip():
            raise ValueError('initial profile requires coordinate/reference-frame provenance')
        if not isinstance(initial_bands,list) or len(initial_bands)!=profile_bands:
            raise ValueError('initial profile band count must match face_profile_bands')
        if any(not isinstance(b,dict) or set(b)!={'z_low_m','z_high_m','front_x_m'} for b in initial_bands):
            raise ValueError('initial bands need exact height and front-coordinate fields')
        face_profile=Profile(g.upstream_edge_x_m,g.ice_back_edge_x_m,g.ice_bottom_z_m,g.crest_z_m,
                             g.extrusion_width_m,tuple(Band(**b) for b in initial_bands))
        if face_profile.floor_z_m!=g.floor_z_m or face_profile.bands[0].front_x_m!=g.face_x_m:
            raise ValueError('initial profile and basal routing geometry disagree')
    elif 'profile_provenance' in upgrades:raise ValueError('profile provenance supplied without initial bands')
    keys={'bypass_floor_m','bypass_floor_limit_m','bypass_width_m','bypass_area_m2','weir_coefficient_sqrt_m_s',
          'obstruction_height_m','aperture_m','discharge_coefficient','darcy_conductance_m2_s',
          'water_h_W_m2_K','face_cover_fraction','debris_capacity_J_m2_K','debris_conductance_W_m2_K',
          'initial_debris_temperature_C','new_debris_temperature_C','initial_surface_cold_J_m2',
          'contact_manning_n','contact_friction_slope',
          'new_surface_cold_J_m2','ice_subcooling_K','ice_heat_capacity_J_kg_K','max_meltwater_fraction'}
    if not isinstance(p,dict) or set(p)!=keys:raise ValueError('unknown/missing physics fields')
    signed={'bypass_floor_m','bypass_floor_limit_m','initial_debris_temperature_C','new_debris_temperature_C'}
    for k,v in p.items():finite(v,k,None if k in signed else 0)
    toe_config=config_from_dict(upgrades.get('toe_heat_transfer'),legacy_fixed_h_W_m2_K=p['water_h_W_m2_K'])
    bypass_config=config_from_dict(upgrades.get('bypass_heat_transfer'),legacy_fixed_h_W_m2_K=p['water_h_W_m2_K'])
    bypass_slope=upgrades.get('bypass_friction_slope')
    if bypass_slope is not None and finite(bypass_slope,'bypass friction slope',0)<=0:raise ValueError('positive bypass friction slope required')
    if bypass_config.method!='fixed' and bypass_slope is None:
        raise ValueError('nonfixed bypass transfer requires its own explicit friction slope')
    for k in ('bypass_width_m','bypass_area_m2','weir_coefficient_sqrt_m_s','debris_capacity_J_m2_K','ice_heat_capacity_J_kg_K','contact_manning_n','contact_friction_slope'):
        if p[k]<=0:raise ValueError(k+' must be positive')
    f=p['face_cover_fraction']
    if not 0<=f<1:raise ValueError('face cover fraction must be >=0 and <1')
    fixed_cover=upgrades.get('fixed_cover_interval_m')
    if fixed_cover is not None:
        if not isinstance(fixed_cover,list) or len(fixed_cover)!=2:raise ValueError('fixed cover requires two elevations')
        lo,hi=[finite(v,'fixed cover elevation') for v in fixed_cover]
        if not g.floor_z_m<=lo<hi<=g.crest_z_m:raise ValueError('fixed cover must lie within the initial face')
        if abs((hi-lo)/(g.crest_z_m-g.floor_z_m)-f)>1e-12:raise ValueError('initial cover fraction and fixed footprint disagree')
        if not isinstance(upgrades.get('cover_interval_provenance'),str) or not upgrades['cover_interval_provenance'].strip():
            raise ValueError('fixed cover needs coordinate/mask provenance')
    elif 'cover_interval_provenance' in upgrades:raise ValueError('cover provenance supplied without fixed mask')
    if wet_mode!='reject' and f>0:
        if any(p[k]!=0 for k in ('ice_subcooling_K','initial_surface_cold_J_m2','new_surface_cold_J_m2')):
            raise ValueError('explicit wet-debris modes require temperate ice and zero extra cold reservoirs')
        if min(p['initial_debris_temperature_C'],p['new_debris_temperature_C'])<0:
            raise ValueError('wet-debris modes do not implement cold debris/refreezing')
        if wet_mode=='resolved_temperate' and p['debris_conductance_W_m2_K']<=0:
            raise ValueError('resolved diffusion requires positive debris conductance')
    if face_profile is not None and any(p[k]!=0 for k in ('ice_subcooling_K','initial_surface_cold_J_m2','new_surface_cold_J_m2')):
        raise ValueError('height-resolved face adapter requires temperate ice and zero cold reservoirs')
    if p['bypass_floor_limit_m']>=p['bypass_floor_m']:raise ValueError('bypass domain invalid')
    dtmax=finite(c['step_s'],'step_s',0)
    if dtmax==0:raise ValueError('positive step required')
    if not isinstance(c['periods'],list) or not c['periods']:raise ValueError('forcing periods required')
    end=0.;periods=[]
    for row in c['periods']:
        if set(row)!={'duration_s','inflow_m3_s','water_superheat_K','face_atmospheric_net_W_m2',
                     'floor_atmospheric_net_W_m2','bypass_atmospheric_net_W_m2'}:raise ValueError('forcing schema')
        for k,v in row.items():finite(v,k,0 if k in ('duration_s','inflow_m3_s','water_superheat_K') else None)
        if face_profile is not None and row['face_atmospheric_net_W_m2']<0:
            raise ValueError('height profile currently requires nonnegative prescribed face heat; no spatial cold-state solver')
        if row['duration_s']==0:raise ValueError('positive forcing duration required')
        periods.append((end,end+row['duration_s'],row));end+=row['duration_s']
    removal=c['removal_time_s']
    if removal is not None:
        finite(removal,'removal time',0)
        if removal>end:raise ValueError('removal beyond horizon')
    if math.ceil(end/dtmax)>200000:raise ValueError('too many steps')
    par=Parameters(p['debris_capacity_J_m2_K'],p['debris_conductance_W_m2_K'],L)
    def cover_bounds(cliff):
        if fixed_cover is not None:return tuple(fixed_cover)
        h=f*(cliff.crest_z_m-cliff.floor_z_m)
        return (cliff.crest_z_m-h,cliff.crest_z_m) if cover_location=='upper' else (cliff.floor_z_m,cliff.floor_z_m+h)
    def covered_area(cliff):
        lo,hi=cover_bounds(cliff)
        return max(0.,min(hi,cliff.crest_z_m)-max(lo,cliff.floor_z_m))*cliff.extrusion_width_m
    initial_cover_area=covered_area(g)
    areas={'face_bare':g.face_area_m2*(1-f) if fixed_cover is None else g.face_area_m2-initial_cover_area,'floor':g.floor_area_m2,'bypass':p['bypass_area_m2']}
    if f>0:areas['face_debris']=g.face_area_m2*f if fixed_cover is None else initial_cover_area
    states={k:State(a,par.capacity_j_m2_k*a*p['initial_debris_temperature_C'] if k=='face_debris' else 0.,
                    a*p['initial_surface_cold_J_m2']) for k,a in areas.items()}
    columns=(ColumnPartition(par.capacity_j_m2_k,0.,states['face_debris'].area_m2,None,
             (p['initial_debris_temperature_C'],)*debris_cells)
             if wet_mode=='resolved_temperate' and f>0 else None)
    def stored():return sum(s.debris_energy_j-s.cold_debt_j for s in states.values())
    initial_stored=stored();initial_ice=face_profile.ice_volume_m3 if face_profile is not None else g.ice_volume_m3;zb=p['bypass_floor_m'];zb0=zb
    totals=dict(atmospheric_J=0.,water_in_J=0.,water_out_J=0.,latent_J=0.,fresh_ice_warming_J=0.,
                remap_enthalpy_J=0.,face_melt_m3=0.,floor_melt_m3=0.,bypass_melt_m3=0.,water_in_m3=0.,
                water_out_m3=0.,phase_water_m3=0.,covered_face_melt_m3=0.,bare_face_melt_m3=0.,
                absolute_surface_forcing_J=0.,absolute_material_enthalpy_transfer_J=0.)
    cold_vol=RHO_I*p['ice_heat_capacity_J_kg_K']*p['ice_subcooling_K'];rows=[];t=0.
    def phase(k,res):
        volume=res.latent_j/(RHO_I*L+cold_vol)
        latent=res.latent_j if cold_vol==0 else volume*RHO_I*L
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
            height=g.crest_z_m-g.floor_z_m
            wet_height=min(depth,height) if qt>1e-14 else 0.
            if depth>height and qt>1e-14:raise ValueError('water overtops the modelled cliff; extra routing is unsupported')
            covered_height=f*height
            cover_lo,cover_hi=cover_bounds(g)
            wet_covered_height=max(0.,min(cover_hi,contact_stage)-max(cover_lo,g.floor_z_m)) if qt>1e-14 else 0.
            wet_debris=wet_covered_height*g.extrusion_width_m
            wet=max(0.,wet_height-wet_covered_height)*g.extrusion_width_m
            area_normalization=None
            if face_profile is not None:
                if profile_alignment=='surface_classes':
                    face_profile=split_at_elevations(face_profile,[contact_stage,cover_lo,cover_hi],max_profile_bands)
                parts=class_partition(face_profile,contact_stage,[(cover_lo,cover_hi)] if f>0 else [])['totals']
                raw_bare=parts['bare_water_m2']+parts['bare_air_m2']
                wet=states['face_bare'].area_m2*(parts['bare_water_m2']/raw_bare) if raw_bare>0 else 0.
                raw_covered=parts['debris_water_m2']+parts['debris_air_m2']
                wet_debris=(states['face_debris'].area_m2*(parts['debris_water_m2']/raw_covered)
                            if f>0 and raw_covered>0 else 0.)
                area_normalization=dict(raw_partition_m2=parts,
                    bare_scale=states['face_bare'].area_m2/raw_bare if raw_bare>0 else None,
                    covered_scale=states['face_debris'].area_m2/raw_covered if f>0 and raw_covered>0 else None)
                if raw_bare==0 or f>0 and raw_covered==0:
                    raise RuntimeError('profile and material-area partition disagree')
                if any(abs(v-1)>1e-12 for v in (area_normalization['bare_scale'],area_normalization['covered_scale']) if v is not None):
                    raise RuntimeError('profile/material area mismatch exceeds roundoff normalization')
            if wet_debris>0 and wet_mode=='reject':
                raise ValueError('submerged face debris requires a separate thermal boundary closure')
            if columns is not None:
                columns=remap_columns(columns,states['face_debris'].area_m2,wet_debris,
                    (p['new_debris_temperature_C'],)*debris_cells).partition
            toe_transfer=local_transfer(qt,g.extrusion_width_m,depth,toe_config)
            # A fixed bypass coefficient does not imply a supplied bypass depth.
            if bypass_slope is None:
                bydepth=None;bytransfer=None;byh=bypass_config.fixed_h_W_m2_K*bypass_config.gain
            else:
                bydepth=contact_depth(qb,p['bypass_width_m'],p['contact_manning_n'],bypass_slope)
                bytransfer=local_transfer(qb,p['bypass_width_m'],bydepth,bypass_config)
                byh=bytransfer['h_W_m2_K']
            hwater=toe_transfer['h_W_m2_K']
            heat=water_heat_rates(qt,Tw,hwater*g.floor_area_m2,hwater*wet)
            byheat=water_heat_rates(qb,Tw,byh*p['bypass_area_m2'],0.)
            at={k:0. for k in states};at['face_bare']=forcing['face_atmospheric_net_W_m2']*(states['face_bare'].area_m2-wet)
            at['floor']=forcing['floor_atmospheric_net_W_m2']*g.floor_area_m2 if qt<=1e-14 else 0.
            at['bypass']=forcing['bypass_atmospheric_net_W_m2']*p['bypass_area_m2'] if qb<=1e-14 else 0.
            if f>0:at['face_debris']=forcing['face_atmospheric_net_W_m2']*(states['face_debris'].area_m2-(wet_debris if heat_policy=='dry_only' else 0.))
            water={'floor':heat['floor_water_heat_W'],'face_bare':heat['face_water_heat_W'],
                   'bypass':byheat['floor_water_heat_W'],'face_debris':0.}
            latent={};vol={};wet_result=None
            if f>0 and wet_mode!='reject':
                if p['ice_subcooling_K']!=0 or states['face_debris'].cold_debt_j!=0 or states['face_bare'].cold_debt_j!=0:
                    raise ValueError('wet debris kernel requires temperate ice and no unresolved face cold reservoir')
                debris_state=states['face_debris'];capacity=par.capacity_j_m2_k*debris_state.area_m2
                # Floor is upstream; bare and covered face patches exchange with
                # the SAME remaining stream, not duplicated inlet-water heat.
                floor_only=water_heat_rates(qt,Tw,hwater*g.floor_area_m2,0.)
                face_inlet=floor_only['outgoing_water_heat_W']/(RHO_W*CP*qt) if qt>0 else 0.
                if columns is not None:
                    wet_result=advance_columns(columns,par.conductance_w_m2_k,dt,face_inlet,RHO_W*CP*qt,
                        hwater*wet,hwater*wet_debris if water_coupling=='coupled' else 0.,forcing['face_atmospheric_net_W_m2'],heat_policy,substep_fraction)
                    columns=wet_result.partition
                else:
                    wet_result=wet_advance(ThermalState(debris_state.debris_energy_j/capacity),
                        ThermalParameters(capacity,par.conductance_w_m2_k*debris_state.area_m2,hwater*wet,hwater*wet_debris if water_coupling=='coupled' else 0.),
                        dt,face_inlet,RHO_W*CP*qt,atmospheric_power_W=at['face_debris'])
                states['face_debris']=replace(debris_state,debris_energy_j=capacity*wet_result.state.temperature_C)
                water['floor']=floor_only['floor_water_heat_W']
                water['face_bare']=wet_result.bare_ice_heat_J/dt
                # Debit the returned interface energy once through the ordinary
                # phase accounting, not as an extra incoming heat source.
                heat=dict(floor_only,outgoing_water_heat_W=wet_result.water_outlet_heat_J/dt)
                latent['face_debris'],vol['face_debris']=phase('face_debris',
                    SimpleNamespace(latent_j=wet_result.covered_ice_heat_J,state=states['face_debris']))
            for k,state in list(states.items()):
                if k=='face_debris' and wet_result is not None:continue
                res=step(state,par,dt,at[k]/state.area_m2,0.) if k=='face_debris' else step_bare(state,dt,at[k]/state.area_m2,water[k]/state.area_m2)
                if abs(res.residual_j)>1e-7*max(1.,abs(res.atmospheric_j)+abs(res.water_j)):raise RuntimeError('port energy failure')
                latent[k],vol[k]=phase(k,res)
            face_latent=latent['face_bare']+latent.get('face_debris',0.)
            if face_profile is None:
                g,geom=advance_latent(g,face_latent,latent['floor'])
            else:
                covered_total=latent.get('face_debris',0.)
                covered_water=(wet_result.wet_covered_ice_heat_J if columns is not None else
                    covered_total*(wet_debris/states['face_debris'].area_m2) if f>0 else 0.)
                covered_dry=(wet_result.dry_covered_ice_heat_J if columns is not None else covered_total-covered_water)
                class_energy=dict(bare_air=at['face_bare']*dt,bare_water=water['face_bare']*dt,
                    covered_dry=covered_dry,covered_water=covered_water)
                if abs(sum(class_energy.values())-face_latent)>1e-8*max(1.,face_latent):
                    raise RuntimeError('face-class latent energy does not match total phase change')
                cover_intervals=[(cover_lo,cover_hi)] if f>0 else []
                face_profile,geom,proxy=advance_with_partition(face_profile,class_energy,latent['floor'],contact_stage,cover_intervals)
                g=SteppedCliff(**proxy).validate()
            zb-=vol['bypass']/p['bypass_area_m2']
            if zb<=p['bypass_floor_limit_m']:raise ValueError('bypass ice domain exhausted')
            phasewater=RHO_I*sum(vol.values())/RHO_W
            ratio=phasewater/(qin*dt) if qin>0 else (None if phasewater>0 else 0.)
            toe_phasewater=RHO_I*(vol['face_bare']+vol.get('face_debris',0.)+vol['floor'])/RHO_W
            bypass_phasewater=RHO_I*vol['bypass']/RHO_W
            toe_ratio=toe_phasewater/(qt*dt) if qt>0 else (None if toe_phasewater>0 else 0.)
            bypass_ratio=bypass_phasewater/(qb*dt) if qb>0 else (None if bypass_phasewater>0 else 0.)
            local_feedback=(toe_ratio is None or bypass_ratio is None or
                (toe_ratio is not None and toe_ratio>p['max_meltwater_fraction']) or
                (bypass_ratio is not None and bypass_ratio>p['max_meltwater_fraction']))
            if local_feedback and melt_policy is None:
                raise ValueError('local route meltwater feedback requires a declared downstream free-drainage contract')
            if melt_policy is not None and phasewater/dt>melt_policy['max_generated_discharge_m3_s']:
                raise ValueError('generated meltwater exceeds declared downstream drainage capacity')
            if qin==0 and phasewater>0:
                if melt_policy is None:raise ValueError('zero-inflow melt requires an explicit downstream free-drainage contract')
                if phasewater/dt>melt_policy['max_generated_discharge_m3_s']:raise ValueError('generated meltwater exceeds declared downstream drainage capacity')
            elif ratio is not None and ratio>p['max_meltwater_fraction']:
                raise ValueError('meltwater hydraulic feedback exceeds declared reduced-model scope')
            newcover=covered_area(g)
            newareas={'face_bare':g.face_area_m2*(1-f) if fixed_cover is None else g.face_area_m2-newcover,'floor':g.floor_area_m2,'bypass':p['bypass_area_m2']}
            if f>0:newareas['face_debris']=g.face_area_m2*f if fixed_cover is None else newcover
            for k,s in list(states.items()):
                if newareas[k]==s.area_m2:continue
                if newareas[k]<s.area_m2:raise RuntimeError('unimplemented material contraction')
                content=AreaContent(par.capacity_j_m2_k*p['new_debris_temperature_C'] if k=='face_debris' else 0.,
                                    p['new_surface_cold_J_m2'],0.)
                rem=remap_area(s,newareas[k],incoming=content);states[k]=rem.state
                totals['remap_enthalpy_J']+=rem.debris_transfer_j-rem.cold_debt_transfer_j
                totals['absolute_material_enthalpy_transfer_J']+=abs(rem.debris_transfer_j-rem.cold_debt_transfer_j)
                if k=='face_debris' and columns is not None:
                    col_remap=remap_columns(columns,newareas[k],columns.wet_area_m2,
                        (p['new_debris_temperature_C'],)*debris_cells)
                    if abs(col_remap.imported_energy_J-rem.debris_transfer_j)>1e-8*max(1.,abs(rem.debris_transfer_j)):
                        raise RuntimeError('debris column/material import mismatch')
                    columns=col_remap.partition
                    if abs(columns.stored_energy_J-states[k].debris_energy_j)>1e-8*max(1.,abs(states[k].debris_energy_j)):
                        raise RuntimeError('debris reservoir/column storage mismatch')
            totals['atmospheric_J']+=sum(at.values())*dt
            totals['absolute_surface_forcing_J']+=sum(abs(v) for v in at.values())*dt
            totals['water_in_J']+=RHO_W*CP*(qt+qb)*Tw*dt
            totals['water_out_J']+=(heat['outgoing_water_heat_W']+byheat['outgoing_water_heat_W'])*dt
            totals['water_in_m3']+=qin*dt
            totals['water_out_m3']+=(qt+qb)*dt+phasewater
            totals['phase_water_m3']+=phasewater
            totals['face_melt_m3']+=vol['face_bare']+vol.get('face_debris',0.)
            totals['covered_face_melt_m3']+=vol.get('face_debris',0.);totals['bare_face_melt_m3']+=vol['face_bare']
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
            rows[-1].update(wet_covered_face_area_m2=wet_debris,face_cover_location=cover_location if fixed_cover is None else 'fixed_elevation_interval',
                covered_z_interval_m=(cover_lo,cover_hi),covered_face_melt_m3=totals['covered_face_melt_m3'],
                toe_heat_transfer=toe_transfer,bypass_heat_transfer=bytransfer,
                wet_debris_thermal=None if wet_result is None else asdict(wet_result),
                covered_heat_policy=heat_policy,debris_thermal_mode=wet_mode,
                water_debris_coupling=water_coupling,
                generated_meltwater_discharge_m3_s=phasewater/dt,
                phase_water_destination='downstream outlet; not recirculated into upstream warm inflow')
            profile_snapshot=(asdict(face_profile) if face_profile is not None and
                              (len(rows)==1 or len(rows)%snapshot_stride==0 or t>=end) else None)
            if face_profile is not None and ledger_detail=='summary':
                detailed={'old_profile_partition','band_class_energy_J','band_face_energy_J','after_ledges','per_band_retreat_m','band_retreat_m','face_band_retreat_m'}
                profile_ledger={k:v for k,v in geom.items() if k not in detailed}
            else:profile_ledger=geom if face_profile is not None else None
            rows[-1].update(face_profile=profile_snapshot,height_profile_active=face_profile is not None,
                surface_partition_sample_time_s=t-dt,state_sample_time_s=t,
                debris_columns_post_geometry=None if columns is None else asdict(columns),
                area_roundoff_normalization=area_normalization,profile_alignment=profile_alignment,
                profile_height_band_count=None if face_profile is None else len(face_profile.bands),
                toe_phase_water_to_route_inflow_ratio=toe_ratio,bypass_phase_water_to_route_inflow_ratio=bypass_ratio,
                local_phase_water_feedback_exceeds_small_source_bound=local_feedback,
                downstream_drainage_contract=None if melt_policy is None else dict(melt_policy),
                cumulative_energy_accounting_magnitude_J=max(1.,abs(initial_stored)+totals['water_in_J']+totals['absolute_surface_forcing_J']+totals['absolute_material_enthalpy_transfer_J']),
                profile_geometry_ledger=profile_ledger,profile_ledger_detail=ledger_detail,
                geometry_is_routing_proxy=face_profile is not None,
                undercut_diagnostic=None if face_profile is None else face_profile.undercut_diagnostic(),
                horizontal_ledge_heat_policy=None if face_profile is None else 'adiabatic, prescribed approximation; no collapse model')
    remaining_ice=face_profile.ice_volume_m3 if face_profile is not None else g.ice_volume_m3
    phase_geom=initial_ice-remaining_ice+(zb0-zb)*p['bypass_area_m2']
    volumes=totals['face_melt_m3']+totals['floor_melt_m3']+totals['bypass_melt_m3']
    scale=max(1.,abs(totals['atmospheric_J'])+totals['water_in_J'])
    if abs(residual)>1e-9*scale or abs(phase_geom-volumes)>1e-9*max(1.,volumes):raise RuntimeError('global conservation failure')
    return dict(schema_version=2,input_sha256=hashlib.sha256(json.dumps(c,sort_keys=True,allow_nan=False).encode()).hexdigest(),
                evidence='conditional coupled reduced-model forecast; not spatial CFD or field validation',
                provenance=c['provenance'],series=rows,totals=totals,
                checks=dict(global_energy_residual_J=residual,geometry_phase_volume_residual_m3=phase_geom-volumes,
                            water_bookkeeping_residual_m3=totals['water_out_m3']-totals['water_in_m3']-totals['phase_water_m3']),
                scope=['thin effective upstream obstacle; bare downstream floor',
                       'explicit upper/lower prescribed face cover; optional normal diffusion with conservative wet/dry pools',
                       ('quasi-steady common stage; height-band face with explicitly adiabatic ledges' if face_profile is not None else 'quasi-steady common stage and uniform planar retreat'),
                       'downstream face wetting from imposed-friction-slope rectangular Manning rating, not upstream stage',
                       'prescribed lumped face fractions; no tracked bare/debris boundary advection',
                       'explicit incoming surface reservoirs and bulk fresh-ice sensible heat',
                       'resolved option uses separate wet/dry normal columns; donor-average class transfer is not full material advection',
                       'meltwater leaves downstream; small-feedback bound enforced',
                       'prescribed net atmospheric face flux; no inferred runoff or water temperature',
                       'no collapse, glacier advection or resolved cold-ice PDE'])

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args()
    if args.input.resolve()==args.output.resolve():parser.error('output cannot replace input')
    try:result=forecast(json.loads(args.input.read_text()))
    except (ValueError,TypeError,KeyError) as exc:
        print(json.dumps(dict(status='REFUSED',reason=str(exc))));raise SystemExit(2)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    import os,tempfile
    with tempfile.NamedTemporaryFile(mode='w',dir=args.output.parent,delete=False) as f:
        temporary=Path(f.name);f.write(json.dumps(result,indent=2,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
    os.replace(temporary,args.output);print(json.dumps(result['checks'],indent=2))
