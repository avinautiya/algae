"""Isolated depth-integrated meltwater source; NOT an integrated flow solver.

Sensible heat is referenced to the melting point. The caller has already
removed latent heat and supplies an ice-volume-equivalent vertical melt depth.
Stationary ice supplies no horizontal momentum and zero sensible enthalpy.
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Result:
    bed_m: np.ndarray
    depth_m: np.ndarray
    momentum_x_m2_s: np.ndarray
    momentum_y_m2_s: np.ndarray
    depth_temperature_m_K: np.ndarray
    ice_volume_lost_m3: float
    water_volume_added_m3: float
    phase_mass_residual_kg: float
    sensible_heat_change_J: float
    maximum_surface_drop_m: float


def apply(*, bed_m, depth_m, momentum_x_m2_s, momentum_y_m2_s,
          depth_temperature_m_K, ice_melt_depth_m, cell_area_m2,
          rho_i_kg_m3=917., rho_w_kg_m3=1000., cp_w_J_kg_K=4217.):
    """Return fresh arrays; never mutate the parent model or reuse its ledgers.

    No drainage occurs in this isolated source step. Rain, refreezing, moving
    ice, porous storage, momentum/pressure work and full ALE balance are outside
    this operator. Active cells must be explicit: set all other melt depths0.
    """
    scalars=[rho_i_kg_m3,rho_w_kg_m3,cp_w_J_kg_K]
    if any(isinstance(x,(bool,np.bool_)) or not np.isscalar(x) or not np.isfinite(x) or x<=0 for x in scalars):
        raise ValueError('densities and heat capacity must be finite positive numbers')
    if rho_i_kg_m3>rho_w_kg_m3:
        raise ValueError('this freshwater-ice operator requires ice density <= water density')
    vals=[np.asarray(x,dtype=float) for x in [bed_m,depth_m,momentum_x_m2_s,momentum_y_m2_s,depth_temperature_m_K,ice_melt_depth_m]]
    z,h,mx,my,ht,dz=vals
    if not z.size or any(x.shape!=z.shape or not np.isfinite(x).all() for x in vals):
        raise ValueError('finite nonempty arrays of identical shape required')
    a=np.asarray(cell_area_m2,dtype=float)
    if a.shape not in [(),z.shape] or not np.isfinite(a).all() or (a<=0).any():
        raise ValueError('cell areas must be finite positive scalar or same-shape array')
    if (h<0).any() or (dz<0).any() or (ht<0).any():
        raise ValueError('negative depth, melt or sensible heat unsupported')
    dry=h==0
    if (mx[dry]!=0).any() or (my[dry]!=0).any() or (ht[dry]!=0).any():
        raise ValueError('dry cells cannot contain momentum or sensible heat')
    liquid_increment=rho_i_kg_m3/rho_w_kg_m3*dz
    new_h=h+liquid_increment;new_z=z-dz
    if not np.isfinite(new_h).all() or not np.isfinite(new_z).all():
        raise ValueError('source update overflows')
    vi=float(np.sum(a*dz));vw=float(np.sum(a*liquid_increment))
    residual=rho_w_kg_m3*vw-rho_i_kg_m3*vi
    if not all(np.isfinite(x) for x in [vi,vw,residual]):
        raise ValueError('phase ledger overflows')
    return Result(new_z.copy(),new_h.copy(),mx.copy(),my.copy(),ht.copy(),vi,vw,float(residual),0.,float(np.max((1-rho_i_kg_m3/rho_w_kg_m3)*dz)))


def route_source_ratio(*, generated_water_m3, interval_s, route_discharge_m3_s):
    """Undefined at zero external route flow; no global-ratio substitution."""
    x=[generated_water_m3,interval_s,route_discharge_m3_s]
    if any(isinstance(v,(bool,np.bool_)) or not np.isscalar(v) or not np.isfinite(v) for v in x):
        raise ValueError('finite scalars required')
    if generated_water_m3<0 or interval_s<=0 or route_discharge_m3_s<0:
        raise ValueError('nonnegative volumes/discharge and positive interval required')
    return None if route_discharge_m3_s==0 else generated_water_m3/(interval_s*route_discharge_m3_s)
