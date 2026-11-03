"""Passive two-layer debris plus cold-ice phase model for DEVELOPMENT studies.

Cell-to-cell/half-cell conductances integrate dz/k(z) exactly, so the material
interface need not coincide with a mesh face. Constant shared debris capacity.
No contact resistance, wet-debris phase enthalpy, moving mesh or ice advection.
"""
import math


def factor(lower,diag,upper):
    b=list(diag);multipliers=[]
    for i in range(1,len(b)):
        x=lower[i-1]/b[i-1];multipliers.append(x);b[i]-=x*upper[i-1]
    return multipliers,b,list(upper)


def solve(fact,rhs):
    factors,b,upper=fact;r=list(rhs)
    for i in range(1,len(r)):r[i]-=factors[i-1]*r[i-1]
    out=[0.]*len(r);out[-1]=r[-1]/b[-1]
    for i in range(len(r)-2,-1,-1):out[i]=(r[i]-upper[i]*out[i+1])/b[i]
    return out


def interpolation(z,y,locations):
    out=[]
    for p in locations:
        if p<z[0] or p>z[-1]:raise ValueError('extrapolation refused')
        for i in range(len(z)-1):
            if z[i]<=p<=z[i+1]:
                a=(p-z[i])/(z[i+1]-z[i]);out.append(y[i]*(1-a)+y[i+1]*a);break
    return out


def predict(*,k_upper,k_lower,capacity,depth,layer_depth,initial_depths,initial_temperatures,
            surface_samples,sample_depths,interval_s=1800.,dt_s=300.,debris_cells=32,
            ice_cells=64,ice_depth=2.,ice_capacity=1925700.,ice_k=2.1):
    for x in [k_upper,k_lower,capacity,depth,layer_depth,interval_s,dt_s,ice_depth,ice_capacity,ice_k]:
        if type(x) not in (int,float) or not math.isfinite(x) or x<=0:raise ValueError('positive finite physical parameters required')
    if not layer_depth<depth:raise ValueError('layer interface must be inside debris')
    if type(debris_cells) is not int or type(ice_cells) is not int or min(debris_cells,ice_cells)<2:raise ValueError('at least2cells per material')
    if not math.isclose(interval_s/dt_s,round(interval_s/dt_s),abs_tol=1e-12):raise ValueError('dt must divide native interval')
    if len(surface_samples)<2 or any(not math.isfinite(x) for x in surface_samples):raise ValueError('complete native boundary required')
    if initial_depths[0]!=0. or initial_depths[-1]!=depth:raise ValueError('initial profile must span exact domain')
    if len(initial_depths)!=len(initial_temperatures) or any(b<=a for a,b in zip(initial_depths,initial_depths[1:])) or any(not math.isfinite(x) for x in initial_temperatures):raise ValueError('valid initial profile required')
    N=debris_cells;M=ice_cells;dz=depth/N;dzi=ice_depth/M;Sd=capacity*dz;Si=ice_capacity*dzi
    centers=[(j+.5)*dz for j in range(N)]
    def conductance(a,b):
        upper=max(0.,min(b,layer_depth)-a) if a<layer_depth else 0.
        upper=min(upper,b-a);return 1/(upper/k_upper+(b-a-upper)/k_lower)
    gd=[conductance(a,b) for a,b in zip(centers,centers[1:])]
    Gtop=conductance(0.,centers[0]);Gd=conductance(centers[-1],depth);Gi=2*ice_k/dzi;gi=ice_k/dzi
    G=1/(1/Gd+1/Gi);rd=Sd/dt_s;ri=Si/dt_s
    dd=[rd+(gd[j-1] if j else Gtop)+(gd[j] if j<N-1 else Gd) for j in range(N)]
    ii=[ri+(gi if j else Gi)+(gi if j<M-1 else Gi) for j in range(M)]
    fd=factor([-x for x in gd],dd,[-x for x in gd]);fi=factor([-gi]*(M-1),ii,[-gi]*(M-1))
    faces=gd+[G]+[gi]*(M-1);rates=[rd]*N+[ri]*M
    cc=[rates[j]+(faces[j-1] if j else Gtop)+(faces[j] if j<N+M-1 else Gi) for j in range(N+M)]
    fc=factor([-x for x in faces],cc,[-x for x in faces])
    Td=interpolation(initial_depths,initial_temperatures,centers);Ti=[0.]*M
    initial_energy=Sd*math.fsum(Td);profiles=[interpolation([0.]+centers+[depth],[surface_samples[0]]+Td+[0.],sample_depths)]
    masses=[];surface_energy=[];deep_energy=[];maximum_defect=0.;cold=0
    steps=int(round(interval_s/dt_s))
    for Ts in surface_samples[:-1]:
        native_mass=[]
        for _ in range(steps):
            rhs=[rd*t for t in Td];rhs[0]+=Gtop*Ts;D=solve(fd,rhs)
            I=solve(fi,[ri*t for t in Ti]);qd=Gd*D[-1];qi=-Gi*I[0]
            if qd>=qi:Tb=0.;latent=(qd-qi)*dt_s
            else:
                rhs=[rd*t for t in Td]+[ri*t for t in Ti];rhs[0]+=Gtop*Ts
                allT=solve(fc,rhs);D,I=allT[:N],allT[N:]
                Tb=(Gd*D[-1]+Gi*I[0])/(Gd+Gi);latent=0.;cold+=1
                if Tb>0.:raise ValueError('positive cold-interface temperature')
            if max(I)>0:raise ValueError('ice above melting')
            Qin=Gtop*(Ts-D[0])*dt_s;Qdeep=Gi*I[-1]*dt_s
            change=Sd*math.fsum(b-a for a,b in zip(Td,D))+Si*math.fsum(b-a for a,b in zip(Ti,I))
            maximum_defect=max(maximum_defect,abs(change+latent-Qin+Qdeep))
            native_mass.append(latent/333700.);surface_energy.append(Qin);deep_energy.append(Qdeep);Td,Ti=D,I
        masses.append(math.fsum(native_mass));profiles.append(interpolation([0.]+centers+[depth],[Ts]+Td+[Tb],sample_depths))
    mass=math.fsum(masses);final_energy=Sd*math.fsum(Td)+Si*math.fsum(Ti)
    residual=final_energy-initial_energy+mass*333700.-math.fsum(surface_energy)+math.fsum(deep_energy)
    return dict(predicted_temperature_C=profiles,melt_mmwe=mass,interval_melt_mmwe=masses,
                energy_residual_J_m2=residual,max_step_energy_residual_J_m2=maximum_defect,
                relative_energy_residual=residual/max(abs(initial_energy),abs(final_energy),mass*333700.,abs(math.fsum(surface_energy)),1.),
                cold_interface_steps=cold,geometry_fixed=True,wet_pore_enthalpy_included=False)
