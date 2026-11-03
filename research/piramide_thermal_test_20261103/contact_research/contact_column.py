"""Conservative fixed-geometry debris/cold-ice conduction with melt complementarity.

Interface temperature Tb<=0 C, melt mass M>=0, Tb*M=0. Negative heat does
not create imaginary refreezing water: it changes sensible cold content in
the explicitly represented ice. This is NOT a moving-front Stefan solver.
"""
import math
if __package__:
    from .linear import tridiagonal as _tridiagonal
else:
    from linear import tridiagonal as _tridiagonal


def num(x,name,positive=False):
    if type(x) not in (int,float) or not math.isfinite(x) or (positive and x<=0):
        raise ValueError(name+' must be finite numeric'+(' and positive' if positive else ''))
    return float(x)


def entropy_time_defect(old,new,storage):
    """Backward-Euler entropy defect C[log(r)-1+1/r]>=0; units J/m²/K."""
    x=(new-old)/(old+273.15)
    if abs(x)<1e-4:
        value=x*x*(.5+x*(-2/3+x*(.75+x*(-.8+x*(5/6)))))
    else:value=math.log1p(x)-x/(1+x)
    return storage*value


def edge_entropy(a,b,g):
    """Conductive edge irreversibility, g(a-b)^2/(absolute a*absolute b)."""
    return g*(a-b)**2/((a+273.15)*(b+273.15))


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
             record_profiles=False, contact_resistance_m2K_W=0.):
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
    R=contact_resistance_m2K_W
    if type(R) not in (int,float) or math.isnan(R) or R<0:
        raise ValueError('contact resistance must be nonnegative finite or positive infinity')
    if type(record_profiles) is not bool:raise ValueError('record_profiles must be bool')
    Td=[num(x,'initial debris temperature') for x in initial_debris_temperature_C]
    Ti=[num(x,'initial ice temperature') for x in initial_ice_temperature_C]
    surface=[num(x,'surface temperature') for x in surface_temperature_C]
    if not Td or not Ti or not surface:raise ValueError('nonempty explicit profiles and forcing required')
    if any(x>0 for x in Ti):raise ValueError('initial ice above melting point unsupported')
    if any(x<=-273.15 for x in Td+Ti+surface):raise ValueError('temperature must exceed absolute zero')
    dzd=dd/len(Td);dzi=di/len(Ti);Sd=Cd*dzd;Si=Ci*dzi
    gd=kd/dzd;gi=ki/dzi;Gtop=2*gd;Gi=2*gi
    Gd=0. if math.isinf(R) else (Gtop if R==0 else 1/(1/Gtop+R))
    G=0. if Gd==0. else 1/(1/Gd+1/Gi)
    if not all(math.isfinite(x) and x>0 for x in [dzd,dzi,Sd,Si,gd,gi,Gtop,Gi]):raise ValueError('unrepresentable grid/material coefficients')
    E0=Sd*math.fsum(Td)+Si*math.fsum(Ti)
    S0=Sd*math.fsum(math.log1p(t/273.15) for t in Td)+Si*math.fsum(math.log1p(t/273.15) for t in Ti)
    entropy_physical=[];entropy_temporal=[];entropy_boundary=[];entropy_step_defects=[]
    surface_steps=[];deep_steps=[];latent_steps=[];step_defects=[];history=[];intervals=[]
    for native_i,Ts in enumerate(surface):
        elapsed=0.;native_melt=[];native_surface=[];native_deep=[]
        while elapsed<cadence:
            h=min(dt,cadence-elapsed)
            if h+elapsed==elapsed:raise ValueError('step cannot advance time')
            # Active melting candidate with common interface temperature0C.
            D=block(Td,Sd,gd,Gtop,Gd,Ts,0.,h)
            I=block(Ti,Si,gi,Gi,Gi,0.,0.,h)
            qd=Gd*D[-1];qi=-Gi*I[0];net=qd-qi
            if net>=0.:
                Tb=0.;mass=net*h/L;branch='melting'
            else:
                # No phase change. Eliminate the zero-capacity interface using
                # its two half-cell thermal resistances in series.
                D,I=coupled(Td,Ti,Sd,Si,gd,gi,Gtop,G,Gi,Ts,h)
                Tb=(Gd*D[-1]+Gi*I[0])/(Gd+Gi)
                qd=Gd*(D[-1]-Tb);qi=Gi*(Tb-I[0])
                mass=0.;branch='cold_interface'
                if Tb>0.:raise ValueError('complementarity active set produced positive cold interface')
            if any(x>0. for x in I):raise ValueError('ice exceeds melting temperature')
            qsurface=Gtop*(Ts-D[0]);qdeep=Gi*I[-1]
            Tface=D[-1]-qd/Gtop
            contact_jump=Tface-Tb
            entropy=0. if R==0 or math.isinf(R) else R*qd*qd/((Tface+273.15)*(Tb+273.15))
            deltaD=Sd*math.fsum(b-a for a,b in zip(Td,D))
            deltaI=Si*math.fsum(b-a for a,b in zip(Ti,I))
            Ein=h*qsurface;Edeep=h*qdeep;latent=L*mass
            defect=deltaD+deltaI+latent-Ein+Edeep
            physical_entropy_rate=math.fsum([edge_entropy(Ts,D[0],Gtop),
                *[edge_entropy(a,b,gd) for a,b in zip(D,D[1:])],
                edge_entropy(D[-1],Tface,Gtop),entropy,edge_entropy(Tb,I[0],Gi),
                *[edge_entropy(a,b,gi) for a,b in zip(I,I[1:])],edge_entropy(I[-1],0.,Gi)])
            temporal=math.fsum([*[entropy_time_defect(a,b,Sd) for a,b in zip(Td,D)],
                                *[entropy_time_defect(a,b,Si) for a,b in zip(Ti,I)]])
            deltaS=Sd*math.fsum(math.log1p((b-a)/(a+273.15)) for a,b in zip(Td,D))+Si*math.fsum(math.log1p((b-a)/(a+273.15)) for a,b in zip(Ti,I))
            boundaryS=Ein/(Ts+273.15)-Edeep/273.15
            entropy_defect=deltaS+latent/273.15-boundaryS-h*physical_entropy_rate-temporal
            entropy_physical.append(h*physical_entropy_rate);entropy_temporal.append(temporal);entropy_boundary.append(boundaryS);entropy_step_defects.append(entropy_defect)
            if not all(math.isfinite(x) for x in [Tb,qd,qi,mass,Ein,Edeep,latent,defect]):raise ValueError('unrepresentable energy/interface result')
            native_surface.append(Ein);native_deep.append(Edeep);native_melt.append(mass)
            surface_steps.append(Ein);deep_steps.append(Edeep);latent_steps.append(latent);step_defects.append(defect)
            elapsed+=h;Td,Ti=D,I
            row=dict(time_s=native_i*cadence+elapsed,step_s=h,branch=branch,
                     interface_temperature_C=Tb,melt_kg_m2=mass,
                     debris_contact_face_temperature_C=Tface,contact_temperature_jump_C=contact_jump,
                     contact_resistance_m2K_W=(R if math.isfinite(R) else 'insulated_limit'),contact_entropy_generation_W_m2K=entropy,
                     contact_law_residual_C=(contact_jump-R*qd if math.isfinite(R) else None),
                     debris_interface_heat_W_m2=qd,ice_interface_heat_W_m2=qi,
                     surface_heat_J_m2=Ein,deep_heat_out_J_m2=Edeep,
                     debris_storage_change_J_m2=deltaD,ice_storage_change_J_m2=deltaI,
                     latent_heat_J_m2=latent,energy_residual_J_m2=defect)
            row.update(physical_entropy_generation_W_m2K=physical_entropy_rate,
                       backward_euler_entropy_defect_J_m2K=temporal,
                       entropy_balance_residual_J_m2K=entropy_defect)
            if record_profiles:row.update(debris_temperature_C=list(Td),ice_temperature_C=list(Ti))
            history.append(row)
        intervals.append(dict(start_s=native_i*cadence,end_s=(native_i+1)*cadence,
                              melt_kg_m2=math.fsum(native_melt),melt_mmwe=math.fsum(native_melt),
                              surface_heat_J_m2=math.fsum(native_surface),deep_heat_out_J_m2=math.fsum(native_deep)))
    Ef=Sd*math.fsum(Td)+Si*math.fsum(Ti)
    latent=math.fsum(latent_steps);supply=math.fsum(surface_steps);deep=math.fsum(deep_steps)
    residual=Ef-E0+latent-supply+deep
    Sf=Sd*math.fsum(math.log1p(t/273.15) for t in Td)+Si*math.fsum(math.log1p(t/273.15) for t in Ti)
    entropy_residual=Sf-S0+latent/273.15-math.fsum(entropy_boundary)-math.fsum(entropy_physical)-math.fsum(entropy_temporal)
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
                latent_heat_J_kg=L,contact_resistance_m2K_W=(R if math.isfinite(R) else 'insulated_limit'),
                contact_entropy_generation_J_m2K=math.fsum(x['step_s']*x['contact_entropy_generation_W_m2K'] for x in history),
                sensible_entropy_change_J_m2K=Sf-S0,latent_entropy_J_m2K=latent/273.15,
                reservoir_entropy_in_J_m2K=math.fsum(entropy_boundary),
                physical_entropy_generation_J_m2K=math.fsum(entropy_physical),
                backward_euler_entropy_defect_J_m2K=math.fsum(entropy_temporal),
                entropy_balance_residual_J_m2K=entropy_residual,
                max_step_entropy_balance_residual_J_m2K=max(abs(x) for x in entropy_step_defects),
                fixed_geometry=True,moving_ice_advection=False,
                wet_debris_phase_change=False,refreezing_water_created=False,
                field_validation_established=False,
                scope='fixed-interface sensible debris, explicit contact resistance and cold-ice complementarity; causal contact interpretation unverified')
