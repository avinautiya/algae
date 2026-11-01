"""Conservative serial floor-to-face water heat allocation at melting point.

An explicit reduced closure, not a replacement for a resolved temperature field.
"""
import math


def water_heat_rates(flow_m3_s,superheat_K,floor_UA_W_K,face_UA_W_K,
                     rho_water=1000.,cp_water=4217.):
    values=locals().copy()
    for name,v in values.items():
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0:
            raise ValueError(name+' must be finite nonnegative numeric')
    if rho_water==0 or cp_water==0:raise ValueError('water properties must be positive')
    capacity=rho_water*cp_water*flow_m3_s
    incoming=capacity*superheat_K
    floor=incoming*(-math.expm1(-floor_UA_W_K/capacity)) if capacity>0 else 0.
    remaining=incoming-floor
    face=remaining*(-math.expm1(-face_UA_W_K/capacity)) if capacity>0 else 0.
    out=remaining-face
    return dict(incoming_water_heat_W=incoming,floor_water_heat_W=floor,
                face_water_heat_W=face,outgoing_water_heat_W=out,
                outgoing_superheat_K=out/capacity if capacity>0 else 0.,
                energy_residual_W=incoming-floor-face-out,
                scope='serial plug-flow heat exchange at Tm; no freezing or thermal boundary layer resolution')
