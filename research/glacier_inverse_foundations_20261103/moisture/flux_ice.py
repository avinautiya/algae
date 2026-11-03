"""Signed prescribed-flux ice conduction and nonnegative phase-change mass.

Optional zero-thickness *soil* areal heat capacity at the interface is a separate
state. No ice half-cell storage is added: the ice FV cells already contain it.
This is a conditional measured-energy retrieval, not a weather forecast.
"""
import math


def number(x,name,positive=False,nonnegative=False):
    if type(x) not in (int,float) or not math.isfinite(x) or (positive and x<=0) or (nonnegative and x<0):
        raise ValueError(name+' must be finite with the declared sign')
    return float(x)


def tridiagonal(lower,diag,upper,rhs):
    b=list(diag);r=list(rhs)
    for i in range(1,len(b)):
        x=lower[i-1]/b[i-1];b[i]-=x*upper[i-1];r[i]-=x*r[i-1]
    result=[0.]*len(b);result[-1]=r[-1]/b[-1]
    for i in range(len(b)-2,-1,-1):result[i]=(r[i]-upper[i]*result[i+1])/b[i]
    if not all(math.isfinite(x) for x in result):raise ValueError('nonfinite solve')
    return result


def soil_storage_components(T4_C,T5_C,*,z4_m,z5_m,ice_depth_m,soil_capacity_J_m3K):
    """Exact P1 sensible storage below the midpoint gradient-flux plane.

    For a=(z5-z4)/2, g=d-z5, Ebelow=Cs[a*T4/4+(3a/4+g/2)*T5]+Cs*g/2*Gamma.
    Returns the observed F and separate interface coefficient A; do not put the
    Gamma term into F again. Reference energy is at0C.
    """
    z4=number(z4_m,'z4');z5=number(z5_m,'z5');d=number(ice_depth_m,'ice depth');C=number(soil_capacity_J_m3K,'soil capacity',positive=True)
    T4=number(T4_C,'T4');T5=number(T5_C,'T5')
    if not z4<z5<=d:raise ValueError('ordered probes and ice below deepest probe required')
    a=(z5-z4)/2;gap=d-z5
    return dict(observed_storage_F_J_m2=C*(a*T4/4+(3*a/4+gap/2)*T5),
                interface_areal_capacity_J_m2K=C*gap/2,
                flux_plane_depth_m=(z4+z5)/2,soil_gap_m=gap,
                piecewise_linear_storage=True,ice_storage_not_included=True)


def prepare_flux(plane_flux_W_m2,observed_storage_F_J_m2,*,interval_seconds):
    """Conservative source conversion with explicit endpoints and no filling.

    qplane is already a declared interval representative. F has one more value
    than intervals; its full difference is removed exactly, not clipped.
    """
    cadence=number(interval_seconds,'interval',positive=True)
    q=[number(x,'plane flux') for x in plane_flux_W_m2]
    F=[number(x,'observed below-plane storage') for x in observed_storage_F_J_m2]
    if not q or len(F)!=len(q)+1:raise ValueError('F endpoints must fully cover every flux interval')
    return dict(effective_interface_flux_W_m2=[x-(b-a)/cadence for x,a,b in zip(q,F,F[1:])],
                original_plane_flux_W_m2=q,observed_storage_F_J_m2=F,
                interval_seconds=cadence,storage_change_J_m2=F[-1]-F[0])


def simulate(*,interface_flux_W_m2,initial_ice_temperature_C,interval_seconds,dt_seconds,
             ice_depth_m=2.,ice_conductivity_W_mK=2.1,ice_capacity_J_m3K=917.*2100.,
             initial_interface_temperature_C=0.,soil_interface_areal_capacity_J_m2K=0.,
             latent_heat_J_kg=333700.,record_profiles=False):
    """Flux-positive-down boundary, ice deep temperature0C, fixed geometry.

    Γ<=0C, M>=0, Γ*M=0. q can have either sign. Cooling stores sensible cold
    content in ice and the optional soil interface capacity; no water refreezes.
    """
    q=[number(x,'signed prescribed flux') for x in interface_flux_W_m2]
    T=[number(x,'initial ice temperature') for x in initial_ice_temperature_C]
    if not q or not T:raise ValueError('nonempty explicit forcing and ice profile required')
    Gamma=number(initial_interface_temperature_C,'initial interface temperature')
    if Gamma>0 or max(T)>0:raise ValueError('ice/interface cannot start above melting')
    if min([Gamma]+T)<=-273.15:raise ValueError('initial ice/interface must exceed absolute zero')
    A=number(soil_interface_areal_capacity_J_m2K,'soil interface areal capacity',nonnegative=True)
    cadence=number(interval_seconds,'interval',positive=True);step=number(dt_seconds,'dt',positive=True)
    depth=number(ice_depth_m,'ice depth',positive=True);k=number(ice_conductivity_W_mK,'ice conductivity',positive=True)
    C=number(ice_capacity_J_m3K,'ice capacity',positive=True);L=number(latent_heat_J_kg,'latent heat',positive=True)
    if type(record_profiles) is not bool:raise ValueError('record_profiles must be bool')
    n=len(T);dz=depth/n;S=C*dz;g=k/dz;G=2*g
    E0=S*math.fsum(T)+A*Gamma;input_energy=[];deep_energy=[];masses=[];defects=[];history=[];intervals=[]
    for j,flux in enumerate(q):
        elapsed=0.;native_mass=[];native_deep=[]
        while elapsed<cadence:
            h=min(step,cadence-elapsed)
            if h+elapsed==elapsed:raise ValueError('dt cannot advance time')
            rate=S/h
            # Active phase candidate: interface at melting point, including
            # energy needed to warm a previously cold soil-interface state.
            diag=[rate+(g if i else G)+(g if i<n-1 else G) for i in range(n)]
            warm=tridiagonal([-g]*(n-1),diag,[-g]*(n-1),[rate*t for t in T])
            qi=-G*warm[0];excess=flux-qi-A*(0.-Gamma)/h
            if excess>=0.:
                updated=warm;newGamma=0.;mass=excess*h/L;branch='melting'
            elif A==0.:
                # Eliminate Γ: qinput=G(Γ-T0) -> Γ=T0+qinput/G.
                diag=[rate+(g if i else 0.)+(g if i<n-1 else G) for i in range(n)]
                rhs=[rate*t for t in T];rhs[0]+=flux
                updated=tridiagonal([-g]*(n-1),diag,[-g]*(n-1),rhs)
                newGamma=updated[0]+flux/G;mass=0.;branch='cold_neumann'
            else:
                # Explicit boundary soil state plus full ice block. The
                # interface unknown carries A, never an extra ice half cell.
                faces=[G]+[g]*(n-1);rates=[A/h]+[rate]*n
                diag=[rates[i]+(faces[i-1] if i else 0.)+(faces[i] if i<n else G) for i in range(n+1)]
                rhs=[A*Gamma/h+flux]+[rate*t for t in T]
                solution=tridiagonal([-x for x in faces],diag,[-x for x in faces],rhs)
                newGamma,updated=solution[0],solution[1:];mass=0.;branch='cold_interface_capacity'
            if newGamma>0. or max(updated)>0.:raise ValueError('phase active set exceeded melting temperature')
            if min([newGamma]+updated)<=-273.15:raise ValueError('forcing drove ice below absolute zero')
            Qin=flux*h;Qdeep=G*updated[-1]*h
            ice_change=S*math.fsum(b-a for a,b in zip(T,updated));sheet_change=A*(newGamma-Gamma)
            latent=L*mass;defect=ice_change+sheet_change+latent-Qin+Qdeep
            if not all(math.isfinite(x) for x in [newGamma,Qin,Qdeep,ice_change,sheet_change,latent,defect]):raise ValueError('unrepresentable phase/energy result')
            input_energy.append(Qin);deep_energy.append(Qdeep);masses.append(mass);defects.append(defect)
            native_mass.append(mass);native_deep.append(Qdeep);elapsed+=h;Gamma,T=newGamma,updated
            row=dict(time_s=j*cadence+elapsed,step_s=h,branch=branch,
                     interface_temperature_C=Gamma,melt_kg_m2=mass,
                     prescribed_flux_W_m2=flux,ice_storage_change_J_m2=ice_change,
                     soil_interface_storage_change_J_m2=sheet_change,deep_heat_out_J_m2=Qdeep,
                     latent_heat_J_m2=latent,energy_residual_J_m2=defect)
            if record_profiles:row['ice_temperature_C']=list(T)
            history.append(row)
        intervals.append(dict(start_s=j*cadence,end_s=(j+1)*cadence,melt_mmwe=math.fsum(native_mass),
                              supplied_heat_J_m2=flux*cadence,deep_heat_out_J_m2=math.fsum(native_deep)))
    Ef=S*math.fsum(T)+A*Gamma;mass=math.fsum(masses);supplied=math.fsum(input_energy);deep=math.fsum(deep_energy)
    residual=Ef-E0+L*mass-supplied+deep
    return dict(melt_mmwe=mass,final_interface_temperature_C=Gamma,final_ice_temperature_C=T,
        initial_sensible_storage_J_m2=E0,final_sensible_storage_J_m2=Ef,
        supplied_heat_J_m2=supplied,deep_heat_out_J_m2=deep,latent_heat_J_m2=L*mass,
        energy_residual_J_m2=residual,relative_energy_residual=residual/max(abs(E0),abs(Ef),abs(supplied),L*mass,abs(deep),1.),
        max_step_energy_residual_J_m2=max(abs(x) for x in defects),
        complementarity_max_product=max(abs(x['interface_temperature_C']*x['melt_kg_m2']) for x in history),
        history=history,intervals=intervals,soil_interface_areal_capacity_J_m2K=A,
        signed_flux_clipped=False,refreezing_water_created=False,geometry_fixed=True,
        ice_advection_included=False,wet_soil_phase_enthalpy_included=False,
        field_accuracy_established=False)


def retrieve_from_plane(*,plane_flux_W_m2,observed_storage_F_J_m2,**ice_options):
    """Full ledger adds observed below-plane soil storage exactly once."""
    forcing=prepare_flux(plane_flux_W_m2,observed_storage_F_J_m2,interval_seconds=ice_options['interval_seconds'])
    result=simulate(interface_flux_W_m2=forcing['effective_interface_flux_W_m2'],**ice_options)
    plane_heat=math.fsum(q*forcing['interval_seconds'] for q in forcing['original_plane_flux_W_m2'])
    changeF=forcing['storage_change_J_m2']
    residual=(result['final_sensible_storage_J_m2']-result['initial_sensible_storage_J_m2']+changeF+
              result['latent_heat_J_m2']-plane_heat+result['deep_heat_out_J_m2'])
    result.update(original_plane_heat_J_m2=plane_heat,observed_soil_storage_change_J_m2=changeF,
                  full_plane_energy_residual_J_m2=residual,
                  full_plane_ledger='Delta(Eice+A*Gamma+F)+L*M=Qplane-Qdeep',
                  observed_F_and_interface_A_counted_once=True)
    return result
