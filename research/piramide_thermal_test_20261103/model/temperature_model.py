"""Native-cadence signed sensible-conduction predictor; no implicit melt claim.

Allows two observed Dirichlet boundaries for process-only validation. Negative
debris temperatures are never clipped. Normalized flux/energy below must NOT
be called melt or joules without an independent bulk heat capacity and an
independently justified ice/phase boundary.
"""
import math
from thermal_memory.solver import _tridiagonal


def number(x, name, positive=False):
    if type(x) not in (int,float) or not math.isfinite(x) or (positive and x<=0):
        raise ValueError(name+' must be finite numeric'+(' and positive' if positive else ''))
    return float(x)


def interpolate(depths, values, locations):
    if len(depths)!=len(values) or len(depths)<2 or any(b<=a for a,b in zip(depths,depths[1:])):
        raise ValueError('ordered depth-value pairs required')
    answer=[]
    for z in locations:
        if z<depths[0] or z>depths[-1]:raise ValueError('profile extrapolation refused')
        for i in range(len(depths)-1):
            if depths[i]<=z<=depths[i+1]:
                f=(z-depths[i])/(depths[i+1]-depths[i])
                answer.append(values[i]*(1-f)+values[i+1]*f);break
    return answer


def predict(alpha_m2_s, *, domain_depth_m, initial_depths_m, initial_temperature_C,
            upper_temperature_C, lower_temperature_C, interval_seconds,
            sample_depths_m, cells, dt_seconds, hold='left'):
    """Return sample temperatures at every native label, including initial time.

    Initial profile includes both actual domain boundaries. `left` uses sample
    i on [ti,ti+1); `right` uses i+1. These are explicit numerical reconstruction
    scenarios of source samples, not assertions about source averaging semantics.
    Both series must span the identical native clock. No gap filling occurs.
    """
    alpha=number(alpha_m2_s,'diffusivity',True)
    depth=number(domain_depth_m,'domain depth',True)
    cadence=number(interval_seconds,'native cadence',True)
    step=number(dt_seconds,'dt',True)
    if type(cells) is not int or cells<2:raise ValueError('integer cells >=2 required')
    if hold not in ('left','right'):raise ValueError('explicit left/right hold required')
    upper=[number(t,'upper temperature') for t in upper_temperature_C]
    lower=[number(t,'lower temperature') for t in lower_temperature_C]
    initial=[number(t,'initial temperature') for t in initial_temperature_C]
    zz=[number(z,'initial depth') for z in initial_depths_m]
    samples=[number(z,'sample depth') for z in sample_depths_m]
    if len(upper)!=len(lower) or len(upper)<2:raise ValueError('same complete native boundary support required')
    if not zz or zz[0]!=0 or zz[-1]!=depth:raise ValueError('initial profile must span exact boundaries')
    if initial[0]!=upper[0] or initial[-1]!=lower[0]:raise ValueError('initial boundary values disagree')
    dz=depth/cells
    centers=[(i+.5)*dz for i in range(cells)]
    T=interpolate(zz,initial,centers)
    def sample(values,top,base):return interpolate([0.]+centers+[depth],[top]+values+[base],samples)
    profiles=[sample(T,upper[0],lower[0])]
    energies=[];basal=[];surface=[];max_defect=0.;subzero=any(t<0 for t in T);negative_basal_substeps=0
    g=alpha/dz;gb=2*g
    initial_store=dz*math.fsum(T)
    for native_i in range(len(upper)-1):
        j=native_i if hold=='left' else native_i+1
        top,base=upper[j],lower[j]
        elapsed=0.;outgoing=[];incoming=[]
        while elapsed<cadence:
            h=min(step,cadence-elapsed)
            if h<=0 or h+elapsed==elapsed:raise ValueError('dt cannot advance native interval')
            rate=dz/h
            diag=[rate+(g if i else gb)+(g if i<cells-1 else gb) for i in range(cells)]
            rhs=[rate*x for x in T];rhs[0]+=gb*top;rhs[-1]+=gb*base
            if not all(math.isfinite(x) for x in diag+rhs):raise ValueError('unrepresentable system')
            updated=_tridiagonal([-g]*(cells-1),diag,[-g]*(cells-1),rhs)
            qi=gb*(top-updated[0]);qo=gb*(updated[-1]-base)
            negative_basal_substeps+=int(qo<0)
            change=dz*math.fsum(b-a for a,b in zip(T,updated))
            defect=change-h*(qi-qo)
            max_defect=max(max_defect,abs(defect));incoming.append(h*qi);outgoing.append(h*qo)
            T=updated;elapsed+=h;subzero=subzero or min(T)<0
        surface.append(math.fsum(incoming));basal.append(math.fsum(outgoing))
        energies.append(dz*math.fsum(T));profiles.append(sample(T,top,base))
    ledger=(energies[-1]-initial_store)-math.fsum(surface)+math.fsum(basal)
    if not all(math.isfinite(x) for x in [ledger,*energies,*surface,*basal]):raise ValueError('unrepresentable ledger')
    return dict(predicted_temperature_C=profiles,
                interval_surface_normalized_heat_K_m=surface,
                interval_basal_normalized_heat_K_m=basal,
                normalized_energy_residual_K_m=ledger,
                maximum_step_normalized_energy_residual_K_m=max_defect,
                native_interval_seconds=cadence,sample_depths_m=samples,
                boundary_reconstruction=hold+'_hold',subzero_debris_computed_without_clipping=subzero,
                negative_basal_flux_substeps=negative_basal_substeps,
                melt_prediction=False,
                scope='signed sensible conduction; pore freezing/phase change not represented')


def conditional_melt(normalized_basal_heat, *, bulk_heat_capacity_J_m3K,
                     latent_heat_J_kg, ice_interface_contract):
    """Explicitly conditional latent conversion; no sensor-plane/ice conflation.

    Contract must be reviewed externally. No negative heat is clipped or silently
    converted to zero melt. Wet-pore phase change and cold ice are unsupported.
    """
    C=number(bulk_heat_capacity_J_m3K,'bulk capacity',True)
    L=number(latent_heat_J_kg,'latent heat',True)
    if not isinstance(ice_interface_contract,dict):raise ValueError('ice-interface contract required')
    for name in ('base_is_ice_interface','basal_temperature_zero_supported',
                 'pore_phase_change_negligible_supported','bulk_capacity_independent_of_test_melt'):
        if ice_interface_contract.get(name) is not True:raise ValueError('unsupported physical contract: '+name)
    if not isinstance(ice_interface_contract.get('source'),str) or not ice_interface_contract['source'].strip():
        raise ValueError('physical contract source required')
    values=[number(x,'normalized basal heat') for x in normalized_basal_heat]
    if not values or any(x<0 for x in values):raise ValueError('negative basal heat/refreezing unsupported; no clipping')
    return dict(conditional_melt_mmwe=C*math.fsum(values)/L,
                conditional_interval_melt_mmwe=[C*x/L for x in values],
                physical_contract=ice_interface_contract,
                field_validation_established=False)
