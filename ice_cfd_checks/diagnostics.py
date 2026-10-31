"""Explicit-SI exported-data diagnostics, with no automatic physics pass/fail."""
from math import isclose, log, sqrt

from .references import DamBreak, rectangular_jump, valid


def norms(errors, weights):
    if not errors or len(errors) != len(weights):
        raise ValueError("Need matched, nonempty errors and positive weights")
    for value in errors:
        valid(value,"error")
    for value in weights:
        valid(value,"weight",0.0,True)
    total=sum(weights)
    return dict(L1=sum(w*abs(e) for e,w in zip(errors,weights))/total,
                L2=sqrt(sum(w*e*e for e,w in zip(errors,weights))/total),
                Linf=max(abs(e) for e in errors))


def score_jump(data):
    """Two remote sections; measured Q at each section, depth/velocity idealization."""
    width=data["width_m"]
    up, down=data["upstream"],data["downstream"]
    gravity=data.get("gravity_m_s2",9.81); rho=data.get("rho_kg_m3",1000.0)
    ref=rectangular_jump(up["depth_m"],up["discharge_m3_s"],width,gravity,rho)
    h2=valid(down["depth_m"],"downstream depth",0.0,True)
    q2=valid(down["discharge_m3_s"],"downstream discharge",0.0,True)
    u2=q2/(width*h2)
    m2=rho*(q2*q2/(width*h2)+0.5*gravity*width*h2*h2)
    measured_loss=ref["specific_energy_up_m"]-h2-u2*u2/(2*gravity)
    return dict(reference=ref, measured_head_loss_m=measured_loss,
        depth_down_error_m=h2-ref["depth_down_m"],
        depth_down_error_relative=(h2-ref["depth_down_m"])/ref["depth_down_m"],
        discharge_difference_m3_s=q2-ref["discharge_m3_s"],
        discharge_difference_relative=q2/ref["discharge_m3_s"]-1,
        ideal_momentum_difference_N=m2-ref["momentum_up_N"],
        ideal_momentum_difference_relative=(m2-ref["momentum_up_N"])/ref["momentum_up_N"],
        head_loss_error_m=measured_loss-ref["head_loss_m"],
        head_loss_error_over_upstream_energy=(measured_loss-ref["head_loss_m"])/ref["specific_energy_up_m"],
        Fr_down=u2/sqrt(gravity*h2),
        limitations=["Assumes steady, horizontal rectangular reach, hydrostatic sections and unit velocity-profile factors",
                     "Use actual integrated momentum/pressure fluxes to diagnose nonuniform sections",
                     "Bed friction and external forces between sections may invalidate the ideal comparison",
                     "This does not validate jump location, turbulence, interface heat transfer or melting"],
        status="DIAGNOSTIC_ONLY_NO_ACCEPTANCE_THRESHOLD")


def score_rest(rows, surface_m, depth_scale_m, gravity=9.81):
    valid(surface_m,"surface")
    valid(depth_scale_m,"reference depth",0.0,True)
    valid(gravity,"gravity",0.0,True)
    errors=[]; speeds=[]; wet_errors=[]
    for row in rows:
        bed=valid(float(row["bed_m"]),"bed")
        h=valid(float(row["depth_m"]),"depth",0.0)
        u=valid(float(row["velocity_m_s"]),"velocity")
        target=max(surface_m-bed,0.0)
        errors.append(h-target); speeds.append(abs(u))
        if target>0:
            wet_errors.append(bed+h-surface_m)
    if not errors:
        raise ValueError("No rest samples")
    return dict(depth_error_m=norms(errors,[1.0]*len(errors)),
        max_wet_surface_error_m=max(map(abs,wet_errors),default=0.0),
        max_speed_m_s=max(speeds), max_speed_over_gravity_scale=max(speeds)/sqrt(gravity*depth_scale_m),
        depth_normalization_m=depth_scale_m,
        note="Samples are equally weighted; this is not a domain-volume conservation test")


def score_dam(rows, reference, time_s):
    valid(time_s,"time",0.0,True)
    eh=[]; eq=[]; widths=[]; previous_right=None
    for row in rows:
        left=float(row["left_m"]); right=float(row["right_m"])
        if previous_right is not None and left < previous_right:
            raise ValueError("Dam-break cells must be sorted and nonoverlapping")
        h=valid(float(row["depth_mean_m"]),"mean depth",0.0)
        q=valid(float(row["hu_mean_m2_s"]),"mean hu")
        exact_h,exact_q=reference.cell_average(left,right,time_s)
        widths.append(right-left); eh.append(h-exact_h); eq.append(q-exact_q)
        previous_right=right
    return dict(depth_error_m=norms(eh,widths), hu_error_m2_s=norms(eq,widths),
        depth_scale_m=reference.left_depth_m,
        hu_scale_m2_s=reference.left_depth_m*sqrt(reference.gravity*reference.left_depth_m),
        time_s=time_s, reference_wave_speeds_m_s=reference.speeds,
        note="Compares cell means; long-wave inviscid hydrostatic reference, not exact near-gate VOF dynamics")


def integrate_faces(rows, rho=1000.0, cp=4217.0):
    """Liquid-only outward fluxes; actual numerical liquid flux is required.

Temperature must be the donor/face value consistent with the solver's enthalpy
flux. A reconstructed alpha*U*area is reported for audit, never substituted.
"""
    valid(rho,"rho",0.0,True); valid(cp,"cp",0.0,True)
    qin=qout=ein=eout=area=weighted_temp=reconstructed=0.0
    count=0
    for r in rows:
        a=valid(float(r["area_m2"]),"area",0.0,True)
        alpha=valid(float(r["liquid_fraction"]),"liquid fraction",0.0)
        if alpha>1:
            raise ValueError("Liquid fraction above 1; do not silently clip VOF errors")
        u=valid(float(r["normal_velocity_m_s"]),"normal velocity")
        q=valid(float(r["liquid_volume_flux_m3_s"]),"numerical liquid flux")
        temp=valid(float(r["temperature_excess_K"]),"temperature excess")
        area+=alpha*a; weighted_temp+=alpha*a*temp
        reconstructed+=alpha*a*u
        if q>=0:
            qout+=q; eout+=rho*cp*q*temp
        else:
            qin-=q; ein-=rho*cp*q*temp
        count+=1
    if not count:
        raise ValueError("No face samples")
    # Sign-separated averages remain meaningful when net Q cancels in backflow.
    return dict(faces=count, numerical_net_outflow_m3_s=qout-qin,
        reconstructed_alpha_u_area_m3_s=reconstructed,
        reconstructed_minus_numerical_m3_s=reconstructed-(qout-qin),
        gross_outflow_m3_s=qout,gross_inflow_m3_s=qin,
        advective_heat_out_W=eout,advective_heat_in_W=ein,net_advective_heat_out_W=eout-ein,
        outflow_weighted_temperature_K=eout/(rho*cp*qout) if qout else None,
        inflow_weighted_temperature_K=ein/(rho*cp*qin) if qin else None,
        wet_area_weighted_temperature_K=weighted_temp/area if area else None,
        liquid_area_m2=area,
        note="No conduction included; align normal signs and donor temperature with actual solver flux")


THERMAL_COLUMNS=("water_sensible_J","advective_in_J","advective_out_J",
                 "conductive_in_J","conductive_out_J","external_J","to_ice_J")


def score_thermal(rows):
    """Cumulative water-energy ledger; all flux columns in J since a common origin.

to_ice includes total water-to-ice heat, not only latent heat. Conductive patch
sets and to_ice must be disjoint, otherwise the same exchange is counted twice.
"""
    rows=list(rows)
    if not rows:
        raise ValueError("Empty thermal ledger")
    clean=[]; previous=None
    for r in rows:
        t=valid(float(r["time_s"]),"time",0.0)
        if previous is not None and t<=previous:
            raise ValueError("Thermal samples must strictly increase in time")
        clean.append(dict(time_s=t,**{k:valid(float(r[k]),k) for k in THERMAL_COLUMNS}))
        previous=t
    first=clean[0]; residuals=[]
    for r in clean:
        d={k:r[k]-first[k] for k in THERMAL_COLUMNS}
        residual=d["water_sensible_J"]+d["advective_out_J"]+d["conductive_out_J"]+d["to_ice_J"]-d["advective_in_J"]-d["conductive_in_J"]-d["external_J"]
        scale=max(abs(first["water_sensible_J"])+abs(r["water_sensible_J"])
                  +sum(abs(d[k]) for k in THERMAL_COLUMNS if k!="water_sensible_J"),1.0)
        residuals.append(dict(time_s=r["time_s"],residual_J=residual,scale_J=scale,relative=residual/scale))
    return dict(samples=residuals,max_abs_residual_J=max(abs(r["residual_J"]) for r in residuals),
        max_abs_relative=max(abs(r["relative"]) for r in residuals),
        note="Conservation diagnostic only; check all terms, absolute errors and shared-flux independence")


def convergence(levels):
    """Conditional three-grid Richardson diagnostic, with equal spacing ratios."""
    if len(levels)!=3:
        raise ValueError("Provide exactly three levels from coarse to fine")
    h=[valid(float(r["spacing_m"]),"spacing",0.0,True) for r in levels]
    f=[valid(float(r["value"]),"value") for r in levels]
    if not h[0]>h[1]>h[2]:
        raise ValueError("Spacings must decrease")
    ratios=[h[0]/h[1],h[1]/h[2]]
    if not isclose(ratios[0],ratios[1],rel_tol=1e-6):
        raise ValueError("This diagnostic requires equal refinement ratios")
    d0,d1=f[0]-f[1],f[1]-f[2]
    base=dict(values=f,spacings_m=h,successive_differences=[d0,d1],refinement_ratio=ratios[0])
    if d0*d1<=0:
        return dict(base,status="NO_MONOTONE_ORDER",observed_order=None,extrapolated=None)
    order=log(abs(d0/d1))/log(ratios[0])
    if order<=0:
        return dict(base,status="NOT_CONVERGING",observed_order=order,extrapolated=None)
    extrap=f[2]+(f[2]-f[1])/(ratios[0]**order-1)
    return dict(base,status="CONDITIONAL_ASYMPTOTIC_ESTIMATE",observed_order=order,
        extrapolated=extrap,estimated_fine_error=abs(extrap-f[2]),
        note="Not a confidence interval or proof of asymptotic behavior; check timestep/solver error separately")
