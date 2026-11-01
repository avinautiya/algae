"""Necessary exposed-cliff budget certificate, not a complete survival model."""
import math

def capacity(width,height,theta,beta,packing):
    if not all(math.isfinite(v) for v in (width,height,theta,beta,packing)) or not (width>0 and height>0 and 0<beta<theta<=math.pi/2 and 0<packing<=1):
        raise ValueError('invalid planar geometry')
    return packing*width*height**2/2*(1/math.tan(beta)-1/math.tan(theta))

def audit(*,width,height,mantle_solid_m,englacial_solid_fraction,stock_capacity,initial_stock,external_input,export_upper,energy,retreat,rho_i=917.,latent=334000.):
    for k,v in locals().copy().items():
        if isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) or v<0:raise ValueError(k)
    if min(width,height,rho_i,latent)<=0 or englacial_solid_fraction>=1 or initial_stock>stock_capacity:raise ValueError('invalid physical state')
    release=width*(mantle_solid_m+englacial_solid_fraction*height)
    required=max(0.,release*retreat+initial_stock+external_input-stock_capacity)
    available=stock_capacity-initial_stock+export_upper-external_input
    debris_ceiling=max(0.,available/release) if release else (math.inf if available>=0 else 0.)
    energy_ceiling=energy/(rho_i*latent*width*height)
    return dict(required_export_m3=required,debris_retreat_ceiling_m=debris_ceiling,
                energy_retreat_ceiling_m=energy_ceiling,
                status='conditionally feasible' if required<=export_upper and retreat<=energy_ceiling else 'incompatible',
                qualification='necessary exposure-budget test; no field validation')
