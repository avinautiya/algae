"""Auditable freshwater regime screening, not a melt-rate closure or validity proof."""
from math import isfinite, sqrt


def check(value,name,minimum=None,strict=False):
    if not isfinite(value) or (minimum is not None and (value<minimum or (strict and value==minimum))):
        raise ValueError(f"Invalid {name}")
    return value


def density(temperature_C,scale="ITS-90"):
    """UNESCO EOS-80 pure-water polynomial at 1 atm, restricted to 0..40 deg C.

IPTS-68 coefficients; near-ambient ITS-90 conversion t68=1.00024*t90 is explicit.
This is an established approximate EOS, not a modern metrological reference,
salinity/pressure model or a thermodynamically complete water-property package.
"""
    check(temperature_C,"temperature",0.0)
    if temperature_C>40 or scale not in ("ITS-90","IPTS-68"):
        raise ValueError("Use 0..40 deg C and ITS-90 or IPTS-68")
    t=temperature_C*(1.00024 if scale=="ITS-90" else 1.0)
    coefficients=(999.842594,6.793952e-2,-9.095290e-3,1.001685e-4,-1.120083e-6,6.536332e-9)
    result=0.0
    for value in reversed(coefficients):
        result=result*t+value
    return result


def expansion(temperature_C,scale="ITS-90"):
    density(temperature_C,scale)
    factor=1.00024 if scale=="ITS-90" else 1.0
    t=temperature_C*factor
    derivative=(6.793952e-2-2*9.095290e-3*t+3*1.001685e-4*t*t
                -4*1.120083e-6*t**3+5*6.536332e-9*t**4)*factor
    return -derivative/density(temperature_C,scale)


def regime_report(data):
    required={"temperature_C","depth_m","velocity_m_s","length_m","width_m","ice_position"}
    allowed=required|{"local_patch","potential_drop_m","external_power_W","properties"}
    if not required.issubset(data) or set(data)-allowed:
        raise ValueError("Missing or unknown regime inputs")
    t=data["temperature_C"]; rho=density(t); rho0=density(0.0)
    h=check(data["depth_m"],"depth",0.0,True)
    u=check(data["velocity_m_s"],"speed",0.0)
    length=check(data["length_m"],"length",0.0,True)
    width=check(data["width_m"],"width",0.0,True)
    position=data["ice_position"]
    if position not in ("bed","ceiling"):
        raise ValueError("Only horizontal bed/ceiling screening is implemented")
    props=dict(cp_J_kg_K=4217.0,k_W_m_K=0.561,nu_m2_s=1.79e-6,
               surface_tension_N_m=0.075,rho_ice_kg_m3=917.0,L_J_kg=334000.0,g_m_s2=9.81)
    updates=data.get("properties",{})
    if set(updates)-set(props):
        raise ValueError("Unknown material property")
    props.update(updates)
    for key,value in props.items():
        check(value,key,0.0,True)
    cp,k,nu,sigma,ri,latent,g=(props[name] for name in
        ("cp_J_kg_K","k_W_m_K","nu_m2_s","surface_tension_N_m","rho_ice_kg_m3","L_J_kg","g_m_s2"))
    alpha=k/(rho*cp)
    delta=rho-rho0
    # z points upward. Bed: cold/light below warm/heavy for 0<T<4 deg C.
    delta_top_minus_bottom=delta if position=="bed" else -delta
    n2=-g*delta_top_minus_bottom/(rho*h)
    q=width*h*u; adv=rho*cp*q*t
    drop=check(data.get("potential_drop_m",0.0),"potential drop",0.0)
    power=check(data.get("external_power_W",0.0),"external power")
    mechanical=rho*g*q*drop
    results=dict(rho_bulk_kg_m3=rho,rho_at_0C_kg_m3=rho0,
        expansion_coefficient_per_K=expansion(t),kappa_m2_s=alpha,
        Re_Uh=u*h/nu,Pe_Uh=u*h/alpha,Pr=nu/alpha,
        Fr=u/sqrt(g*h),We_h=rho*u*u*h/sigma,Bo_h=rho*g*h*h/sigma,
        Ste_sensible_over_latent=cp*t/latent,buoyancy_N2_per_s2=n2,
        stratification="unstable_mean_density_difference" if n2<0 else "stable_mean_density_difference" if n2>0 else "neutral_mean_density_difference",
        signed_Ri_N2_h2_over_U2=n2*h*h/(u*u) if u else None,
        signed_Ra_destabilizing=-n2*h**4/(nu*alpha),
        transit_time_s=length/u if u else None,
        thermal_diffusion_h2_over_kappa_s=h*h/alpha,
        buoyancy_time_s=1/sqrt(abs(n2)) if n2 else None,
        discharge_m3_s=q,advected_sensible_power_W=adv,available_potential_power_W=mechanical,
        potential_to_advected_power_ratio=mechanical/adv if adv else None,
        external_to_advected_power_ratio=power/adv if adv else None,
        listed_power_latent_ice_volume_equivalent_m3_s=max(adv+mechanical+power,0.0)/(ri*latent))
    warnings=["Screening ratios are not a linear-stability theorem or a validated inclusion/exclusion threshold",
        "Bulk density contrast can miss a nonmonotone interior profile, especially across 4 deg C",
        "Properties are declared constants; use temperature/uncertainty sensitivity before environmental extrapolation",
        "The latent volume equivalent is a power conversion, not a universal upper bound; storage, kinetic/pressure work and other sources need the full ledger",
        "Zero bulk speed makes a bulk Richardson number undefined, not zero"]
    local=data.get("local_patch")
    if local is not None:
        if set(local)!={"depth_m","speed_m_s","length_m"}:
            raise ValueError("Local patch needs depth_m, speed_m_s and length_m")
        hl=check(local["depth_m"],"local depth",0.0,True)
        ul=check(local["speed_m_s"],"local speed",0.0)
        ll=check(local["length_m"],"local length",0.0,True)
        nl=-g*delta_top_minus_bottom/(rho*hl)
        results["local_patch"]=dict(signed_Ri_N2_h2_over_U2=nl*hl*hl/(ul*ul) if ul else None,
            residence_length_over_speed_s=ll/ul if ul else None,
            buoyancy_time_s=1/sqrt(abs(nl)) if nl else None,
            density_assumption="Same endpoint temperatures as bulk; replace with measured/resolved local thermal contrast")
        warnings.append("Local speed/length is a sensitivity input, not a prediction of recirculation residence time")
    return dict(input=data,properties=props,diagnostics=results,warnings=warnings,
                status="SCREENING_ONLY_NO_MELT_FORECAST_OR_NOVELTY_CLAIM")
