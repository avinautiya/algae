"""Pure ice-only lowering with no unreported floor clipping.

This fixes a bed-geometry operator, not a complete flow/phase-change scheme.
Water/momentum sources, reservoir boundaries and multirate storage remain the
coupled solver's responsibility. No automatic monkey-patching is performed.
"""
import math
import numpy as np


def lower_ice(bed, ice, heat_flux_W_m2, dt_s, cell_area_m2, *,
              rho_i=917., latent_J_kg=334000., minimum_ice_elevation_m=None):
    z=np.asarray(bed,dtype=float);mask=np.asarray(ice)
    heat=np.asarray(heat_flux_W_m2,dtype=float)
    if z.ndim!=2 or mask.dtype!=np.bool_ or mask.shape!=z.shape or heat.shape!=z.shape:
        raise ValueError('matching 2-D bed, Boolean ice and heat arrays required')
    if not np.isfinite(z).all() or not np.isfinite(heat).all() or np.any(heat[mask]<0):
        raise ValueError('finite inputs and nonnegative admitted ice heat required')
    vals=(dt_s,cell_area_m2,rho_i,latent_J_kg)
    if not all(not isinstance(v,(bool,np.bool_)) and math.isfinite(v) for v in vals) or dt_s<0 or min(cell_area_m2,rho_i,latent_J_kg)<=0:
        raise ValueError('invalid physical inputs')
    with np.errstate(over='raise',invalid='raise'):
        try:
            dz=np.where(mask,heat*dt_s/(rho_i*latent_J_kg),0.)
            updated=z-dz
        except FloatingPointError as exc:raise ValueError('numeric range exceeded') from exc
    if minimum_ice_elevation_m is not None:
        if isinstance(minimum_ice_elevation_m,bool) or not math.isfinite(minimum_ice_elevation_m):
            raise ValueError('finite declared ice-domain lower bound required')
        if np.any(updated[mask]<minimum_ice_elevation_m):
            raise RuntimeError('ice-domain floor reached; refusing unaccounted heat clipping')
    actual=np.where(mask,z-updated,0.)
    requested_volume=float(dz.sum()*cell_area_m2)
    actual_volume=float(actual.sum()*cell_area_m2)
    return dict(bed=updated,lowering_m=actual,ice_volume_m3=actual_volume,
                heat_supplied_J=float((heat[mask]*dt_s).sum()*cell_area_m2),
                geometric_latent_J=actual_volume*rho_i*latent_J_kg,
                requested_volume_m3=requested_volume,
                nonice_unchanged=bool(np.array_equal(updated[~mask],z[~mask])))
