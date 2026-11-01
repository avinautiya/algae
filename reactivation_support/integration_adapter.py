"""Opt-in replacement for Sim.lower_bed; never installed automatically.

Use ONLY when the model's continuity equation has not already included meltwater
and the surrounding boundary/energy integration supplies the latent heat. This
adds the phase-volume source, not a second latent heat sink. A full cumulative
physical-time energy ledger remains the caller's responsibility.
"""
import numpy as np
from .phase_update import update


def conservative_lower_bed(sim,dtm,*,rho_i=917.,rho_w=1000.,cp=4217.,latent=334000.):
    if isinstance(dtm,bool) or not np.isfinite(dtm) or dtm<=0:raise ValueError('positive finite morphology step required')
    if any(isinstance(v,bool) or not np.isfinite(v) or v<=0 for v in (rho_i,rho_w,cp,latent)):
        raise ValueError('positive finite material constants required')
    dz=sim.last_qacc*dtm/(rho_i*latent)*sim.ice
    r=update(sim.zb,sim.h,sim.hu,sim.hv,sim.hT,dz,liquid_energy_delta_J_m2=0.,rho_i=rho_i,rho_w=rho_w,cp=cp,latent=latent)
    sim.zb[:]=r['zb'];sim.h[:]=r['h'];sim.hu[:]=r['hu'];sim.hv[:]=r['hv'];sim.hT[:]=r['hT']
    sim.phase_water_added_m3=getattr(sim,'phase_water_added_m3',0.)+float(r['water_added_depth'].sum()*sim.dx**2)
    sim.phase_ice_removed_kg=getattr(sim,'phase_ice_removed_kg',0.)+float(dz.sum()*sim.dx**2*rho_i)
    return float(dz.sum()*sim.dx**2),float((dz*sim.toe).sum()*sim.dx**2),float((dz*sim.byp).sum()*sim.dx**2)
