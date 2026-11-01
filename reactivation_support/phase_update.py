"""Conservative liquid-mass source for a melting shallow-water ice floor.

This is a model repair, not a new physical law. Temperatures use Tm=0, hT=h*T.
The caller MUST declare the liquid sensible-energy change. Use zero only when
latent exchange is already supplied/accounted by the surrounding flow solve.
For a closed adiabatic box use delta_E_liquid=-rho_i*Lf*lowering.
No refreezing or moving-solid horizontal momentum is inferred automatically.
"""
import math
import numpy as np


def update(zb,h,hu,hv,hT,lowering,*,liquid_energy_delta_J_m2,
           rho_i=917.,rho_w=1000.,cp=4217.,latent=334000.,melt_u=0.,melt_v=0.):
    constants=(rho_i,rho_w,cp,latent)
    if any(isinstance(v,bool) or not math.isfinite(v) or v<=0 for v in constants):
        raise ValueError('positive finite material constants required')
    arrays=np.broadcast_arrays(*[np.asarray(v,dtype=float) for v in
        (zb,h,hu,hv,hT,lowering,liquid_energy_delta_J_m2,melt_u,melt_v)])
    if any(not np.all(np.isfinite(a)) for a in arrays):raise ValueError('finite fields required')
    z,h0,mx,my,sensible,dz,dE,um,vm=arrays
    if np.any(h0<0) or np.any(sensible<0) or np.any(dz<0):raise ValueError('negative depth, heat content or melt')
    if np.any((h0==0)&((sensible!=0)|(mx!=0)|(my!=0))):raise ValueError('dry cells cannot store liquid heat or momentum')
    water_added=rho_i/rho_w*dz
    new_sensible=sensible+dE/(rho_w*cp)
    if np.any(new_sensible<0):raise ValueError('declared melt/energy change overdraws available sensible heat')
    hn=h0+water_added
    T=np.divide(new_sensible,hn,out=np.zeros_like(hn),where=hn>0)
    if np.any((hn==0)&(new_sensible!=0)):raise ValueError('sensible heat without liquid')
    return dict(zb=(z-dz).copy(),h=hn.copy(),hu=(mx+water_added*um).copy(),
                hv=(my+water_added*vm).copy(),hT=new_sensible.copy(),temperature=T,
                water_added_depth=water_added.copy(),ice_removed_depth=dz.copy(),
                latent_J_m2=(rho_i*latent*dz).copy(),
                liquid_sensible_delta_J_m2=dE.copy(),
                required_external_energy_J_m2=(rho_i*latent*dz+dE).copy(),
                surface_displacement=((rho_i/rho_w-1)*dz).copy())
