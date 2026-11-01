"""A sufficient NO-CONTACT bound for a stated, uniformly lowering toe model.

Not a classifier of possible recovery, a spatial CFD bound, or a novelty claim.
"""
from dataclasses import dataclass, asdict
import math


def finite(value, name, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name+' must be finite numeric')
    if nonnegative and value < 0:raise ValueError(name+' must be nonnegative')
    return float(value)


@dataclass(frozen=True)
class BoundInputs:
    initial_stage_m: float
    bypass_control_elevation_upper_m: float
    bypass_rating_coefficient_lower_m12_s: float
    total_positive_reservoir_inflow_upper_m3_s: float
    initial_uniform_toe_floor_m: float
    minimum_contact_offset_above_ice_m: float
    toe_floor_area_m2: float
    toe_floor_positive_heat_upper_J: float
    toe_floor_cold_content_lower_J: float = 0.
    rho_ice_kg_m3: float = 917.
    latent_heat_J_kg: float = 334000.


def no_contact_bound(inputs):
    """Check Q_b >= k_min [eta-b_max]_+**(3/2) for the WHOLE horizon.

    Required assumptions: fixed positive reservoir storage derivative; all
    positive water sources are bounded by Imax, no negative bypass flux;
    bypass elevation never above b_max and k never below k_min; uniform toe
    lowering; the stated positive-heat bound includes every energy source;
    no collapse, uplift, conduit below the floor, or unrepresented bypass.
    """
    if not isinstance(inputs, BoundInputs):raise TypeError('BoundInputs required')
    d=asdict(inputs)
    signed={'initial_stage_m','bypass_control_elevation_upper_m','initial_uniform_toe_floor_m'}
    for name,value in d.items():finite(value,name,nonnegative=name not in signed)
    for name in ('bypass_rating_coefficient_lower_m12_s','toe_floor_area_m2',
                 'rho_ice_kg_m3','latent_heat_J_kg'):
        if d[name]<=0:raise ValueError(name+' must be positive')
    eq=d['bypass_control_elevation_upper_m']+(d['total_positive_reservoir_inflow_upper_m3_s']/
        d['bypass_rating_coefficient_lower_m12_s'])**(2/3)
    stage=max(d['initial_stage_m'],eq)
    latent=max(0.,d['toe_floor_positive_heat_upper_J']-d['toe_floor_cold_content_lower_J'])
    melt=latent/(d['rho_ice_kg_m3']*d['latent_heat_J_kg']*d['toe_floor_area_m2'])
    contact=d['initial_uniform_toe_floor_m']-melt+d['minimum_contact_offset_above_ice_m']
    margin=contact-stage
    return dict(stage_upper_m=stage,uniform_floor_lowering_upper_m=melt,
                contact_elevation_lower_m=contact,no_contact_margin_m=margin,
                verdict='contact_excluded_under_stated_bounds' if margin>0 else 'undetermined',
                evidence='analytic comparison bound for uniform-toe reduced model; not a field result',
                inputs=d)
