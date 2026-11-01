"""Fully wetted smooth narrow-slot resistance, with entrance head loss.

This separate closure avoids applying an inviscid thin-opening law to a long
microscopic under-rock gap. It is not a turbulent/rough-granular flow model.
"""
import math


def filled_slot_discharge(*,upstream_stage_m,downstream_stage_m,floor_z_m,
                          aperture_m,width_m,length_m,viscosity_pa_s,
                          rho_water_kg_m3,entrance_discharge_coefficient):
    values=locals().copy()
    for k,v in values.items():
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):
            raise ValueError(k+' must be finite numeric')
    for k in ('aperture_m','width_m','length_m','viscosity_pa_s',
              'rho_water_kg_m3','entrance_discharge_coefficient'):
        if values[k]<=0:raise ValueError(k+' must be positive')
    if aperture_m/width_m>.1:raise ValueError('narrow-slot approximation requires aperture/width <=0.1')
    roof=floor_z_m+aperture_m
    if min(upstream_stage_m,downstream_stage_m)<roof:
        raise ValueError('both ends must be fully wetted above the slot roof')
    head=max(0.,upstream_stage_m-downstream_stage_m)
    g=9.80665
    # DeltaH = a Q + b Q² from Poiseuille slot friction + entrance loss.
    a=12*viscosity_pa_s*length_m/(rho_water_kg_m3*g*width_m*aperture_m**3)
    b=1/(2*g*(entrance_discharge_coefficient*width_m*aperture_m)**2)
    q=2*head/(a+math.sqrt(a*a+4*b*head)) if head>0 else 0.
    speed=q/(width_m*aperture_m)
    dh=2*width_m*aperture_m/(width_m+aperture_m)
    return dict(flow_m3_s=q,mean_speed_m_s=speed,
                reynolds_Dh=rho_water_kg_m3*speed*dh/viscosity_pa_s,
                viscous_head_loss_m=a*q,entrance_head_loss_m=b*q*q,
                total_head_m=head,head_residual_m=head-a*q-b*q*q,
                dissipation_W=rho_water_kg_m3*g*q*head,
                scope='fully wetted smooth laminar narrow slot plus assumed entrance loss; check Reynolds/development/roughness separately',
                warning='Dissipation is not automatically latent ice heat; no partially filled, reverse, rough or turbulent solution')
