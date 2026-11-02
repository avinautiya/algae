"""Steady finite-water/solid-resistance reference; not a full cliff model."""
import math

def positive(x,name,zero=False):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<0 or (not zero and x==0):raise ValueError('invalid '+name)
    return float(x)

def prediction(*,water_capacity_W_K,film_conductance_W_K,ice_conductance_W_K,inlet_superheat_K,absorbed_external_heat_W=0.):
    F=positive(water_capacity_W_K,'water capacity');G=positive(film_conductance_W_K,'film conductance')
    K=positive(ice_conductance_W_K,'ice conductance');T=positive(inlet_superheat_K,'superheat',True)
    qa=positive(absorbed_external_heat_W,'absorbed external heat',True)
    H=F*(-math.expm1(-G/F));ratio=1/(1+H/K)
    bare=qa+H*T;covered=bare*ratio
    debris_temperature=(qa+H*T)/(K+H)
    water_exchange=H*(T-debris_temperature)
    return dict(effective_water_conductance_W_K=H,NTU=G/F,film_to_solid_conductance_ratio=G/K,
                covered_ice_power_W=covered,bare_ice_power_W=bare,dry_ice_power_W=qa,
                debris_temperature_C=debris_temperature,water_to_debris_heat_W=water_exchange,
                outgoing_water_temperature_C=T-water_exchange/F,
                wet_minus_dry_ice_power_W=covered-qa,
                cooling_criterion_margin_W=qa-K*T,
                energy_residual_W=qa+water_exchange-covered,
                thermal_reconnection_fraction=ratio,bare_proxy_relative_overprediction=(bare-covered)/covered if covered>0 else None,
                hydraulic_connection='supplied positive flow; not predicted')

def thickness_limit(*,wet_area_m2,conductivity_W_m_K,water_capacity_W_K,film_conductance_W_K,required_fraction):
    A=positive(wet_area_m2,'wet area');k=positive(conductivity_W_m_K,'conductivity')
    F=positive(water_capacity_W_K,'water capacity');G=positive(film_conductance_W_K,'film conductance')
    r=positive(required_fraction,'required fraction')
    if r>1:raise ValueError('fraction exceeds1')
    H=F*(-math.expm1(-G/F))
    return k*A*(1-r)/(r*H)
