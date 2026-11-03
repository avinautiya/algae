"""Conservative fixed-geometry debris/cold-ice conduction with melt complementarity.

Interface temperature Tb<=0 C, melt mass M>=0, Tb*M=0. Negative heat does
not create imaginary refreezing water: it changes sensible cold content in
the explicitly represented ice. This is NOT a moving-front Stefan solver.
"""
import math
from thermal_memory.solver import _tridiagonal


def num(x,name,positive=False):
    if type(x) not in (int,float) or not math.isfinite(x) or (positive and x<=0):
        raise ValueError(name+' must be finite numeric'+(' and positive' if positive else ''))
    return float(x)


def block(T,Cdz,g,gtop,gbottom,top,bottom,h):
    n=len(T);rate=Cdz/h
    diagonal=[rate+(g if i else gtop)+(g if i<n-1 else gbottom) for i in range(n)]
    rhs=[rate*x for x in T];rhs[0]+=gtop*top;rhs[-1]+=gbottom*bottom
    if not all(math.isfinite(x) for x in diagonal+rhs):raise ValueError('unrepresentable system')
    return _tridiagonal([-g]*(n-1),diagonal,[-g]*(n-1),rhs)


def coupled(Td,Ti,Sd,Si,gd,gi,Gtop,Ginterface,Gdeep,Ts,h):
    n=len(Td);m=len(Ti);rates=[Sd/h]*n+[Si/h]*m
    faces=[gd]*(n-1)+[Ginterface]+[gi]*(m-1)
    diagonal=[rates[j]+(faces[j-1] if j else Gtop)+(faces[j] if j<n+m-1 else Gdeep) for j in range(n+m)]
    rhs=[rate*t for rate,t in zip(rates,Td+Ti)];rhs[0]+=Gtop*Ts
    if not all(math.isfinite(x) for x in diagonal+rhs):raise ValueError('unrepresentable coupled system')
    T=_tridiagonal([-g for g in faces],diagonal,[-g for g in faces],rhs)
    return T[:n],T[n:]


def simulate(*,debris_conductivity_W_mK,debris_heat_capacity_J_m3K,debris_depth_m,
             initial_debris_temperature_C,initial_ice_temperature_C,
             surface_temperature_C,interval_seconds,dt_seconds,
             ice_depth_m=2.,ice_conductivity_W_mK=2.1,
             ice_heat_capacity_J_m3K=917.*2100.,latent_heat_J_kg=333700.,
             record_profiles=False):
    """Surface values apply to consecutive complete intervals, left held.

    Temperatures are cell centered. Surface and deep boundaries are Dirichlet;
    the deep ice boundary is fixed at0C. Initial ice must be<=0C, but debris
    and its surface are signed. No temperatures or fluxes are zero clipped.
    Fluxes and melt are per unit common horizontal interface area.
    """
    kd=num(debris_conductivity_W_mK,'debris conductivity',True)
    Cd=num(debris_heat_capacity_J_m3K,'bulk debris capacity',True)
    dd=num(debris_depth_m,'debris thickness',True)
    ki=num(ice_conductivity_W_mK,'ice conductivity',True)
    Ci=num(ice_heat_capacity_J_m3K,'ice capacity',True)
    di=num(ice_depth_m,'ice domain depth',True)
    L=num(latent_heat_J_kg,'latent heat',True)
    cadence=num(interval_seconds,'forcing interval',True);dt=num(dt_seconds,'dt',True)
    if type(record_profiles) is not bool:raise ValueError('record_profiles must be bool')
    Td=[num(x,'initial debris temperature') for x in initial_debris_temperature_C]
    Ti=[num(x,'initial ice temperature') for x in initial_ice_temperature_C]
    surface=[num(x,'surface temperature') for x in surface_temperature_C]
    if not Td or not Ti or not surface:raise ValueError('nonempty explicit profiles and forcing required')
    if any(x>0 for x in Ti):raise ValueError('initial ice above melting point unsupported')
    dzd=dd/len(Td);dzi=di/len(Ti);Sd=Cd*dzd;Si=Ci*dzi
    gd=kd/dzd;gi=ki/dzi;Gd=2*gd;Gi=2*gi;G=1/(1/Gd+1/Gi)
    if not all(math.isfinite(x) and x>0 for x in [dzd,dzi,Sd,Si,gd,gi,Gd,Gi,G]):raise ValueError('unrepresentable grid/material coefficients')
    E0=Sd*math.fsum(Td)+Si*math.fsum(Ti)
    surface_steps=[];deep_steps=[];latent_steps=[];step_defects=[];history=[];intervals=[]
    for native_i,Ts in enumerate(surface):
        elapsed=0.;native_melt=[];native_surface=[];native_deep=[]
        while elapsed<cadence:
            h=min(dt,cadence-elapsed)
            if h+elapsed==elapsed:raise ValueError('step cannot advance time')
            # Active melting candidate with common interface temperature0C.
            D=block(Td,Sd,gd,Gd,Gd,Ts,0.,h)
            I=block(Ti,Si,gi,Gi,Gi,0.,0.,h)
            qd=Gd*D[-1];qi=-Gi*I[0];net=qd-qi
            if net>=0.:
                Tb=0.;mass=net*h/L;branch='melting'
            else:
                # No phase change. Eliminate the zero-capacity interface using
                # its two half-cell thermal resistances in series.
                D,I=coupled(Td,Ti,Sd,Si,gd,gi,Gd,G,Gi,Ts,h)
                Tb=(Gd*D[-1]+Gi*I[0])/(Gd+Gi)
                qd=Gd*(D[-1]-Tb);qi=Gi*(Tb-I[0])
                mass=0.;branch='cold_interface'
                if Tb>0.:raise ValueError('complementarity active set produced positive cold interface')
            if any(x>0. for x in I):raise ValueError('ice exceeds melting temperature')
            qsurface=Gd*(Ts-D[0]);qdeep=Gi*I[-1]
            deltaD=Sd*math.fsum(b-a for a,b in zip(Td,D))
            deltaI=Si*math.fsum(b-a for a,b in zip(Ti,I))
            Ein=h*qsurface;Edeep=h*qdeep;latent=L*mass
            defect=deltaD+deltaI+latent-Ein+Edeep
            if not all(math.isfinite(x) for x in [Tb,qd,qi,mass,Ein,Edeep,latent,defect]):raise ValueError('unrepresentable energy/interface result')
            native_surface.append(Ein);native_deep.append(Edeep);native_melt.append(mass)
            surface_steps.append(Ein);deep_steps.append(Edeep);latent_steps.append(latent);step_defects.append(defect)
            elapsed+=h;Td,Ti=D,I
            row=dict(time_s=native_i*cadence+elapsed,step_s=h,branch=branch,
                     interface_temperature_C=Tb,melt_kg_m2=mass,
                     debris_interface_heat_W_m2=qd,ice_interface_heat_W_m2=qi,
                     surface_heat_J_m2=Ein,deep_heat_out_J_m2=Edeep,
                     debris_storage_change_J_m2=deltaD,ice_storage_change_J_m2=deltaI,
                     latent_heat_J_m2=latent,energy_residual_J_m2=defect)
            if record_profiles:row.update(debris_temperature_C=list(Td),ice_temperature_C=list(Ti))
            history.append(row)
        intervals.append(dict(start_s=native_i*cadence,end_s=(native_i+1)*cadence,
                              melt_kg_m2=math.fsum(native_melt),melt_mmwe=math.fsum(native_melt),
                              surface_heat_J_m2=math.fsum(native_surface),deep_heat_out_J_m2=math.fsum(native_deep)))
    Ef=Sd*math.fsum(Td)+Si*math.fsum(Ti)
    latent=math.fsum(latent_steps);supply=math.fsum(surface_steps);deep=math.fsum(deep_steps)
    residual=Ef-E0+latent-supply+deep
    scale=max(abs(E0),abs(Ef),abs(latent),abs(supply),abs(deep),1.)
    return dict(melt_mmwe=latent/L,melt_kg_m2=latent/L,
                final_debris_temperature_C=Td,final_ice_temperature_C=Ti,
                debris_cell_depths_m=[(j+.5)*dzd for j in range(len(Td))],
                ice_cell_depths_below_interface_m=[(j+.5)*dzi for j in range(len(Ti))],
                initial_sensible_storage_J_m2=E0,final_sensible_storage_J_m2=Ef,
                surface_heat_J_m2=supply,deep_heat_out_J_m2=deep,latent_heat_J_m2=latent,
                energy_residual_J_m2=residual,relative_energy_residual=residual/scale,
                max_step_energy_residual_J_m2=max(abs(x) for x in step_defects),
                complementarity_max_product_C_kg_m2=max(abs(x['interface_temperature_C']*x['melt_kg_m2']) for x in history),
                history=history,intervals=intervals,
                ice_properties=dict(conductivity_W_mK=ki,volumetric_heat_capacity_J_m3K=Ci,depth_m=di,deep_temperature_C=0.),
                debris_properties=dict(conductivity_W_mK=kd,volumetric_heat_capacity_J_m3K=Cd,depth_m=dd),
                latent_heat_J_kg=L,fixed_geometry=True,moving_ice_advection=False,
                wet_debris_phase_change=False,refreezing_water_created=False,
                field_validation_established=False,
                scope='fixed-interface dry/sensible-debris plus cold-ice conduction and complementarity melting')
