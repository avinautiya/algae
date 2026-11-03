"""Conservative signed conduction diagnostics, NOT a validated melt predictor.

Constant k,C and fixed debris thickness d; z=0 at melting ice, z=d at
surface. Temperatures relative to melting. C*T_t=k*T_zz, T(0)=0.
I=int_0^d (1-z/d)*T dz. Then q_ice=(k/d)*Ts-C*dI/dt.
This exact weighted-moment identity is elementary, not an originality claim.
"""
import math
SIGMA=5.670374419e-8

def positive(value,name):
    if isinstance(value,bool) or not math.isfinite(value) or value<=0: raise ValueError(name+' must be finite and positive')
    return value

def radiometric_temperature_C(lwout_W_m2,lwin_W_m2,emissivity,*,surface_is_bare_debris):
    if surface_is_bare_debris is not True: raise ValueError('Observed surface must be independently identified as bare debris; unknown snow is not bare.')
    for x in [lwout_W_m2,lwin_W_m2]:
        if not math.isfinite(x) or x<0: raise ValueError('radiation must be finite and nonnegative')
    if not math.isfinite(emissivity) or not 0<emissivity<=1: raise ValueError('emissivity must lie in (0,1]')
    emitted=lwout_W_m2-(1-emissivity)*lwin_W_m2
    positive(emitted,'emitted radiation')
    return (emitted/(emissivity*SIGMA))**.25-273.15

def conduction_energy_bounds(k_W_mK,d_m,C_J_m3K,surface_temperature_integral_Ks,Tmin_C,Tmax_C):
    for x,n in [(k_W_mK,'k'),(d_m,'d'),(C_J_m3K,'C')]: positive(x,n)
    if not all(math.isfinite(x) for x in [surface_temperature_integral_Ks,Tmin_C,Tmax_C]) or not Tmin_C<=0<=Tmax_C:
        raise ValueError('invalid temperature or integral bounds')
    qs=k_W_mK/d_m*surface_temperature_integral_Ks
    storage=C_J_m3K*d_m*(Tmax_C-Tmin_C)/2
    return {'signed_energy_lower_J_m2':qs-storage,'signed_energy_upper_J_m2':qs+storage,
            'quasisteady_signed_energy_J_m2':qs,'maximum_storage_correction_J_m2':storage,
            'melt_prediction':False,'assumptions':'fixed d, constant k and C, basal T=0, no internal heat/advection, temperature bounded throughout'}

def signed_energy_to_net_melt_mmwe(energy_J_m2,*,temperate_ice_confirmed,all_other_interface_terms_zero_confirmed,latent_J_kg=334000):
    if temperate_ice_confirmed is not True or all_other_interface_terms_zero_confirmed is not True:
        raise ValueError('Requires verified temperate interface and complete interface energy balance')
    positive(latent_J_kg,'latent heat')
    if not math.isfinite(energy_J_m2) or energy_J_m2<0:
        raise ValueError('Negative/net-refreezing energy is outside melt-only scope')
    return energy_J_m2/latent_J_kg # kg/m2 equals mm water equivalent, NOT metres of ice
