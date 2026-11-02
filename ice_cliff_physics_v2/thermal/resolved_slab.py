"""Independent cell-centred 1-D debris diffusion reference.

Implicit Euler in time, second-order conductive face differences in space.
Top temperature is algebraically coupled to a finite-heat-capacity stream and
parallel bare-ice contact; bottom ice is fixed at 0 C. No lumped ODE imported.
"""
from dataclasses import dataclass
import math


def _finite(x,name,positive=False):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or (positive and x<=0):
        raise ValueError('invalid '+name)
    return float(x)


@dataclass(frozen=True)
class SlabParameters:
    thickness_m: float
    conductivity_W_m_K: float
    volumetric_capacity_J_m3_K: float
    area_m2: float

    def validate(self):
        for k,v in self.__dict__.items():_finite(v,k,True)
        return self

    @property
    def diffusion_time_s(self):
        return self.volumetric_capacity_J_m3_K*self.thickness_m**2/self.conductivity_W_m_K


@dataclass(frozen=True)
class ArealSlabParameters:
    """Normal-coordinate slab specified without inventing a physical thickness.

    C_A=rhoc*d [J/m2/K], K_A=k/d [W/m2/K]. In xi=x/d, these
    two coefficients fully determine temperature and heat evolution. A physical
    d, k and rhoc cannot be uniquely recovered and are not claimed here.
    """
    capacity_J_m2_K: float
    conductance_W_m2_K: float
    area_m2: float

    def validate(self):
        for k,v in self.__dict__.items():_finite(v,k,True)
        return self

    @property
    def diffusion_time_s(self):
        return self.capacity_J_m2_K/self.conductance_W_m2_K


@dataclass(frozen=True)
class SlabState:
    temperatures_C: tuple

    def validate(self):
        if not self.temperatures_C:raise ValueError('nonempty cell temperatures required')
        for v in self.temperatures_C:
            _finite(v,'cell temperature')
            if v<0:raise ValueError('cold debris/refreezing unsupported')
        return self


def _water_coefficients(flow,gb,gd):
    if flow==0 or gb+gd==0:
        return 0.,0.,0.,0.,0.,1.,0.
    ratio=(gb+gd)/flow
    if ratio<1e-5:
        ph=1-ratio/2+ratio**2/6-ratio**3/24+ratio**4/120
        omp=ratio/2-ratio**2/6+ratio**3/24-ratio**4/120
    else:
        ph=-math.expm1(-ratio)/ratio;omp=1-ph
    alpha=gd/(gb+gd)
    return gd*ph,gd*(gb/(gb+gd)+ph*alpha),gb*ph,gb*omp*alpha,alpha,math.exp(-ratio),-math.expm1(-ratio)


def _tridiagonal(lower,diag,upper,rhs):
    n=len(diag);c=list(upper);b=list(rhs);d=list(diag)
    if not all(math.isfinite(x) for x in lower+diag+upper+rhs):raise ValueError('derived matrix overflow')
    for i in range(1,n):
        if d[i-1]<=0 or not math.isfinite(d[i-1]):raise ArithmeticError('nonpositive/nonfinite diffusion pivot')
        factor=lower[i-1]/d[i-1]
        d[i]-=factor*c[i-1];b[i]-=factor*b[i-1]
    if d[-1]<=0 or not math.isfinite(d[-1]):raise ArithmeticError('nonpositive/nonfinite diffusion pivot')
    out=[0.]*n;out[-1]=b[-1]/d[-1]
    for i in range(n-2,-1,-1):out[i]=(b[i]-c[i]*out[i+1])/d[i]
    return out


def advance(state,parameters,dt_s,inlet_temperature_C,water_capacity_W_K,
            bare_water_conductance_W_K=0.,debris_water_conductance_W_K=0.,
            atmospheric_power_W=0.):
    """One implicit step; all reported fluxes use the same end-step solution.

    Fluxes are extensive W/J. This is a resolved debris slab but not resolved
    flow, moving ice, a cold-ice PDE, phase change in debris or field validation.
    End-step sampling makes the time discretisation first order.
    """
    if not isinstance(state,SlabState) or not isinstance(parameters,(SlabParameters,ArealSlabParameters)):raise TypeError('typed slab state and parameters required')
    state.validate();parameters.validate();dt=_finite(dt_s,'duration',True)
    tw=_finite(inlet_temperature_C,'inlet temperature');flow=_finite(water_capacity_W_K,'water capacity')
    gb=_finite(bare_water_conductance_W_K,'bare conductance');gd=_finite(debris_water_conductance_W_K,'debris conductance');air=_finite(atmospheric_power_W,'atmospheric power')
    if min(tw,flow,gb,gd)<0:raise ValueError('nonnegative water temperature, capacity and conductances required')
    n=len(state.temperatures_C)
    if isinstance(parameters,SlabParameters):
        dx=parameters.thickness_m/n
        cell_capacity=parameters.volumetric_capacity_J_m3_K*parameters.area_m2*dx
        link=parameters.conductivity_W_m_K*parameters.area_m2/dx
    else:
        dx=None
        cell_capacity=parameters.capacity_J_m2_K*parameters.area_m2/n
        link=parameters.conductance_W_m2_K*parameters.area_m2*n
    half=2*link
    if not all(math.isfinite(x) and x>0 for x in [cell_capacity,link,half]):raise ValueError('derived slab coefficients overflow/underflow')
    a,b,bare_in,bare_debris,alpha,decay,omdecay=_water_coefficients(flow,gb,gd)
    source=air+a*tw
    fraction=half/(half+b);surface_loss=fraction*b;surface_source=fraction*source
    diag=[cell_capacity/dt]*n;rhs=[cell_capacity/dt*t for t in state.temperatures_C]
    lower=[-link]*(n-1);upper=[-link]*(n-1)
    for i in range(n):
        diag[i]+=link*(int(i>0)+int(i<n-1))
    diag[0]+=half;diag[-1]+=surface_loss;rhs[-1]+=surface_source
    temp=_tridiagonal(lower,diag,upper,rhs)
    surface=(source+half*temp[-1])/(half+b)
    if not all(math.isfinite(x) for x in temp+[surface]):raise ValueError('nonfinite diffusion solution')
    if min(temp+[surface])<0:raise ValueError('cold debris/refreezing trajectory unsupported')
    incoming_slab=surface_source-surface_loss*temp[-1]
    water_debris=incoming_slab-air
    bare_heat=(bare_in*tw+bare_debris*surface)*dt
    covered_heat=half*temp[0]*dt
    storage=cell_capacity*math.fsum(new-old for new,old in zip(temp,state.temperatures_C))
    water_in=flow*tw*dt
    water_out=flow*(decay*tw+omdecay*alpha*surface)*dt
    # Check the thermal stream expression separately; the top flux below is the
    # consistent discrete conductive boundary flux used by the linear system.
    independent_qd=(a*tw-b*surface)*dt
    advective_exchange=flow*omdecay*(tw-alpha*surface)*dt
    water_formula_residual=(water_debris*dt+bare_heat)-advective_exchange
    energy_residual=air*dt+water_debris*dt-storage-covered_heat
    initial_storage=cell_capacity*math.fsum(state.temperatures_C)
    passivity_margin=initial_storage+water_in+max(air,0.)*dt-bare_heat-covered_heat
    if not all(math.isfinite(x) for x in [incoming_slab,water_debris,bare_heat,covered_heat,storage,water_in,water_out,independent_qd,advective_exchange,water_formula_residual,energy_residual,passivity_margin]):raise ValueError('derived heat inventory overflow')
    return dict(state=SlabState(tuple(temp)),surface_temperature_C=surface,
        atmospheric_heat_J=air*dt,water_inlet_heat_J=water_in,water_outlet_heat_J=water_out,
        water_to_debris_heat_J=water_debris*dt,bare_ice_heat_J=bare_heat,
        covered_ice_heat_J=covered_heat,total_ice_heat_J=bare_heat+covered_heat,
        debris_storage_change_J=storage,global_energy_residual_J=energy_residual,
        top_water_formula_residual_J=water_debris*dt-independent_qd,
        water_exchange_residual_J=water_formula_residual,passivity_margin_J=passivity_margin,
        cell_capacity_J_K=cell_capacity,cell_width_m=dx,cell_width_normalized=1/n,time_discretization='implicit Euler, end-step fluxes')
