"""Allocate a coupled model's net latent class energies to a height profile.

This adapter does not compute heat-transfer physics or invent a horizontal-ledge
boundary condition. The parent must explicitly impose its ledge policy and gate
cold/refreezing cases before calling this temperate, nonnegative-energy adapter.
"""
from dataclasses import asdict
import math
try:
    from .face_profile import Profile, advance, face_partition, finite, RHO_ICE, LATENT
except ImportError:
    from face_profile import Profile, advance, face_partition, finite, RHO_ICE, LATENT

GEOMETRY_KEYS = frozenset(['upstream_edge_x_m','ice_back_edge_x_m','ice_bottom_z_m',
                         'face_x_m','floor_z_m','crest_z_m','extrusion_width_m'])
CLASSES = {
    'bare_air': 'bare_air_m2',
    'bare_water': 'bare_water_m2',
    'covered_dry': 'debris_air_m2',
    'covered_water': 'debris_water_m2',
}


def initialize_profile(geometry, n_bands=32):
    if not isinstance(geometry, dict) or set(geometry) != GEOMETRY_KEYS:
        raise ValueError('exact seven-field SteppedCliff geometry required')
    return Profile.uniform(**geometry, n_bands=n_bands)


def routing_proxy(profile):
    if not isinstance(profile, Profile):
        raise ValueError('validated Profile required')
    return dict(upstream_edge_x_m=profile.upstream_edge_x_m,
                ice_back_edge_x_m=profile.ice_back_edge_x_m,
                ice_bottom_z_m=profile.ice_bottom_z_m,
                face_x_m=profile.bands[0].front_x_m,
                floor_z_m=profile.floor_z_m,
                crest_z_m=profile.crest_z_m,
                extrusion_width_m=profile.extrusion_width_m)


def _allocate(energy, weights):
    """Nonnegative sequential allocation with a final remainder, no energy clip."""
    used=[(i,w) for i,w in enumerate(weights) if w>0]
    if not used:
        if energy>0:
            raise ValueError('positive class energy supplied for zero old-profile class area')
        return [0.]*len(weights)
    suffix=[0.]*len(used)
    for j in range(len(used)-1,-1,-1):
        suffix[j]=math.fsum([used[j][1],suffix[j+1] if j+1<len(used) else 0.])
    allocation=[0.]*len(weights);remaining=energy
    for j,(i,weight) in enumerate(used):
        value=remaining if j==len(used)-1 else remaining*(weight/suffix[j])
        if value<0 or value>remaining:
            raise ArithmeticError('class allocation failed nonnegative energy contract')
        allocation[i]=value;remaining-=value
    if remaining!=0:
        raise ArithmeticError('unallocated class energy')
    return allocation


def advance_with_partition(profile, energies_J, floor_latent_J, water_stage,
                           cover_z_intervals=(), *, rho_ice=RHO_ICE, latent_heat=LATENT):
    """Advance an existing profile with snapshot vertical-class latent energies.

    Energies are NET LATENT, already excluding storage and fresh-ice warming.
    All four class keys must be supplied, including zeros for inactive classes.
    The supplied stage/masks are old-step Eulerian inputs, not solved wetting.
    Each class has a uniform net latent flux within its old exposed area; the
    resulting per-band average recession converges under height refinement.
    """
    if not isinstance(profile, Profile):
        raise ValueError('validated Profile required')
    if not isinstance(energies_J, dict) or set(energies_J)!=set(CLASSES):
        raise ValueError('exact four class energies required: '+','.join(CLASSES))
    energy={k:finite(v,k+' latent energy',nonnegative=True) for k,v in energies_J.items()}
    floor=finite(floor_latent_J,'floor latent energy',nonnegative=True)
    partition=face_partition(profile,water_stage,cover_z_intervals)
    contributions={}
    for name,area_name in CLASSES.items():
        weights=[row[area_name] for row in partition['bands']]
        try:
            contributions[name]=_allocate(energy[name],weights)
        except ValueError as exc:
            raise ValueError(name+': '+str(exc)) from exc
    band_energy=[math.fsum(contributions[name][i] for name in CLASSES)
                 for i in range(len(profile.bands))]
    after,ledger=advance(profile,band_energy,floor,rho_ice,latent_heat)
    allocated=math.fsum(band_energy);requested=math.fsum(energy.values())
    proxy=routing_proxy(after)
    proxy_volume=after.extrusion_width_m*((after.ice_back_edge_x_m-after.upstream_edge_x_m)*
                    (after.floor_z_m-after.ice_bottom_z_m)+(after.ice_back_edge_x_m-proxy['face_x_m'])*
                    (after.crest_z_m-after.floor_z_m))
    ledger.update(class_energy_J=dict(energy),class_areas_before_m2={name:partition['totals'][area] for name,area in CLASSES.items()},
                  old_profile_partition=partition,band_class_energy_J=contributions,
                  band_face_energy_J=band_energy,allocated_face_energy_J=allocated,
                  allocation_residual_J=allocated-requested,
                  before_profile_ice_volume_m3=profile.ice_volume_m3,
                  after_profile_ice_volume_m3=after.ice_volume_m3,
                  routing_proxy_is_physical_volume=False,
                  routing_proxy_ice_volume_m3=proxy_volume,
                  physical_minus_proxy_remaining_volume_m3=after.ice_volume_m3-proxy_volume,
                  ledge_policy='must be explicitly supplied/gated by parent; no ledge energy added here',
                  state_assumptions='temperate ice, no cold debt/refreezing in this adapter; Eulerian mask; supplied hydraulic stage')
    return after,ledger,proxy
