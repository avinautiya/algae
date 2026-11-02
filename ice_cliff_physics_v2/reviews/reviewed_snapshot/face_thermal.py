"""Two conservative wet/dry pools, each resolving normal debris diffusion.

No instantaneous heat mixing across the full face. Within each wet/dry class,
the normal profile is representative; this is not a height-resolved material
advection solver or resolved free-surface flow.
"""
from dataclasses import dataclass, replace
import math
from thermal.column_partition import ColumnPartition
from thermal.resolved_slab import SlabState,ArealSlabParameters,advance as slab_advance
from thermal.wet_debris import ThermalState,ThermalParameters,advance as bare_advance


@dataclass(frozen=True)
class FaceResult:
    state: ThermalState
    partition: ColumnPartition
    bare_ice_heat_J: float
    covered_ice_heat_J: float
    wet_covered_ice_heat_J: float
    dry_covered_ice_heat_J: float
    atmospheric_heat_J: float
    water_inlet_heat_J: float
    water_outlet_heat_J: float
    water_to_debris_heat_J: float
    debris_storage_change_J: float
    global_energy_residual_J: float
    water_exchange_residual_J: float
    thermal_substeps: int


def advance_columns(partition,conductance_W_m2_K,dt_s,Tin,water_capacity_W_K,
                    bare_G_W_K,debris_G_W_K,net_external_W_m2,heat_policy,
                    max_diffusion_time_fraction=.01):
    partition.validate()
    values=[conductance_W_m2_K,dt_s,max_diffusion_time_fraction]
    if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0 for x in values):
        raise ValueError('positive finite slab conductance, duration and substep fraction required')
    if heat_policy not in ('dry_only','prescribed_all_cover'):
        raise ValueError('explicit covered-surface heat policy required')
    if not math.isfinite(net_external_W_m2):raise ValueError('finite net external flux required')
    if max_diffusion_time_fraction>.1:raise ValueError('thermal substep fraction exceeds declared accuracy policy')
    tau=partition.capacity_J_m2_K/conductance_W_m2_K
    if not math.isfinite(tau) or tau<=0:raise ValueError('unrepresentable thermal diffusion time')
    count=dt_s/(tau*max_diffusion_time_fraction)
    if not math.isfinite(count) or count>10000:raise ValueError('too many thermal substeps')
    n=max(1,math.ceil(count))
    dt=dt_s/n;current=partition;initial=partition.stored_energy_J
    totals=dict(bare=0.,wet_ice=0.,dry_ice=0.,air=0.,water_in=0.,water_out=0.,water_debris=0.)
    for _ in range(n):
        if current.wet_area_m2>0:
            p=ArealSlabParameters(current.capacity_J_m2_K,conductance_W_m2_K,current.wet_area_m2)
            qa=net_external_W_m2*current.wet_area_m2 if heat_policy=='prescribed_all_cover' else 0.
            w=slab_advance(SlabState(current.wet_profile_C),p,dt,Tin,water_capacity_W_K,
                bare_water_conductance_W_K=bare_G_W_K,debris_water_conductance_W_K=debris_G_W_K,
                atmospheric_power_W=qa)
            wet_profile=w['state'].temperatures_C
        else:
            # Bare water shares the stream once; a zero-area debris class does
            # not require inventing a slab or an empty-cell energy reservoir.
            w0=bare_advance(ThermalState(0.),ThermalParameters(0.,0.,bare_G_W_K,0.),dt,Tin,water_capacity_W_K)
            w=dict(bare_ice_heat_J=w0.bare_ice_heat_J,covered_ice_heat_J=0.,
                   atmospheric_heat_J=0.,water_inlet_heat_J=w0.water_inlet_heat_J,
                   water_outlet_heat_J=w0.water_outlet_heat_J,water_to_debris_heat_J=0.)
            wet_profile=None
        if current.dry_area_m2>0:
            p=ArealSlabParameters(current.capacity_J_m2_K,conductance_W_m2_K,current.dry_area_m2)
            d=slab_advance(SlabState(current.dry_profile_C),p,dt,0.,0.,
                          atmospheric_power_W=net_external_W_m2*current.dry_area_m2)
            dry_profile=d['state'].temperatures_C
            totals['dry_ice']+=d['covered_ice_heat_J'];totals['air']+=d['atmospheric_heat_J']
        else:dry_profile=None
        totals['bare']+=w['bare_ice_heat_J'];totals['wet_ice']+=w['covered_ice_heat_J']
        totals['air']+=w['atmospheric_heat_J'];totals['water_in']+=w['water_inlet_heat_J']
        totals['water_out']+=w['water_outlet_heat_J'];totals['water_debris']+=w['water_to_debris_heat_J']
        current=replace(current,wet_profile_C=wet_profile,dry_profile_C=dry_profile)
    storage=current.stored_energy_J-initial
    ice=totals['bare']+totals['wet_ice']+totals['dry_ice']
    residual=totals['air']+totals['water_in']-totals['water_out']-storage-ice
    exchange=totals['water_in']-totals['water_out']-totals['water_debris']-totals['bare']
    mean=current.stored_energy_J/(current.capacity_J_m2_K*current.total_area_m2)
    return FaceResult(ThermalState(mean),current,totals['bare'],totals['wet_ice']+totals['dry_ice'],
        totals['wet_ice'],totals['dry_ice'],totals['air'],totals['water_in'],totals['water_out'],
        totals['water_debris'],storage,residual,exchange,n)
