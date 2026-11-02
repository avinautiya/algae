"""Independent fixed-geometry thermal-storage identity and screening bounds.

This is a linear passive-network identity, not a historical-priority claim.
It does not predict route recovery or validate the production implementation.
"""
from dataclasses import dataclass
import math


def finite(x, name, positive=False):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0 or (positive and x <= 0):
        raise ValueError('invalid '+name)
    return float(x)


@dataclass(frozen=True)
class Network:
    capacity_J_K: float
    ice_conductance_W_K: float
    bare_water_conductance_W_K: float
    debris_water_conductance_W_K: float
    water_capacity_W_K: float

    def __post_init__(self):
        for name, x in self.__dict__.items(): finite(x, name, name == 'capacity_J_K')
        if self.ice_conductance_W_K + self.debris_water_conductance_W_K <= 0:
            raise ValueError('no debris relaxation sink')
        if self.water_capacity_W_K == 0 and self.bare_water_conductance_W_K + self.debris_water_conductance_W_K > 0:
            raise ValueError('water conductance needs positive water capacity')

    def coefficients(self):
        gb=self.bare_water_conductance_W_K;gd=self.debris_water_conductance_W_K;g=gb+gd
        if g == 0: return dict(a=0., b=0., d=0., bare_inlet=0., phi=1., NTU=0.)
        x=g/self.water_capacity_W_K
        phi=1. if x == 0 else -math.expm1(-x)/x
        alpha=gd/g
        return dict(a=gd*phi, b=gd*(1-(1-phi)*alpha), d=gb*(1-phi)*alpha,
                    bare_inlet=gb*phi,phi=phi,NTU=x)

    def diagnostics(self):
        c=self.coefficients();n=self.ice_conductance_W_K+c['b']
        if n <= 0: raise ValueError('no relaxation')
        return dict(relaxation_time_s=self.capacity_J_K/n,
                    terminal_storage_to_ice_factor=(self.ice_conductance_W_K+c['d'])/n,
                    NTU=c['NTU'])


def exact_segment(net, temperature_C, inlet_temperature_C, atmosphere_W, duration_s):
    for name,x in [('temperature',temperature_C),('inlet_temperature',inlet_temperature_C),('atmosphere',atmosphere_W),('duration',duration_s)]:finite(x,name,name=='duration')
    c=net.coefficients();n=net.ice_conductance_W_K+c['b'];p=atmosphere_W+c['a']*inlet_temperature_C
    tau=net.capacity_J_K/n;tss=p/n;x=duration_s/tau;one=-math.expm1(-x)
    t1=temperature_C+(tss-temperature_C)*one
    ramp=(x/2-x*x/6+x**3/24-x**4/120+x**5/720) if x<1e-4 else 1-one/x
    integ=temperature_C*duration_s+(tss-temperature_C)*duration_s*ramp
    ice=c['bare_inlet']*inlet_temperature_C*duration_s+(net.ice_conductance_W_K+c['d'])*integ
    qs=c['bare_inlet']*inlet_temperature_C*duration_s+(net.ice_conductance_W_K+c['d'])*tss*duration_s
    water=(c['a']+c['bare_inlet'])*inlet_temperature_C*duration_s-c['a']*integ
    du=net.capacity_J_K*(t1-temperature_C)
    gamma=net.diagnostics()['terminal_storage_to_ice_factor']
    return dict(temperature_C=t1,ice_J=ice,instantaneous_ice_J=qs,water_J=water,storage_change_J=du,
                identity_residual_J=(qs-ice)-gamma*du,
                energy_residual_J=water+atmosphere_W*duration_s-ice-du)


def rk4_segment(net, temperature_C, inlet_temperature_C, atmosphere_W, duration_s, nsteps):
    """Independent explicit ODE integration of T, ice and water heat."""
    if isinstance(nsteps,bool) or not isinstance(nsteps,int) or nsteps<1:raise ValueError('positive integer steps')
    c=net.coefficients();k=net.ice_conductance_W_K;C=net.capacity_J_K
    def rhs(y):
        t=y[0]
        qd=c['a']*inlet_temperature_C-c['b']*t
        qb=c['bare_inlet']*inlet_temperature_C+c['d']*t
        return ((atmosphere_W+qd-k*t)/C, qb+k*t, qb+qd)
    y=(temperature_C,0.,0.);h=duration_s/nsteps
    def plus(a,b,s):return tuple(x+s*z for x,z in zip(a,b))
    for _ in range(nsteps):
        a=rhs(y);b=rhs(plus(y,a,h/2));cc=rhs(plus(y,b,h/2));d=rhs(plus(y,cc,h))
        y=tuple(v+h*(aa+2*bb+2*ccc+dd)/6 for v,aa,bb,ccc,dd in zip(y,a,b,cc,d))
    return dict(temperature_C=y[0],ice_J=y[1],water_J=y[2])


def screening_bound(net, debris_temperature_excursion_K, ice_area_m2, rho_i=917., latent_J_kg=334000.):
    """Maximum fixed-coefficient phase-energy difference from terminal storage.

    Input excursion is an assumed bound, not inferred from data. The area must
    match the footprint whose equivalent lowering is reported.
    """
    finite(debris_temperature_excursion_K,'temperature excursion')
    finite(ice_area_m2,'ice area',True);finite(rho_i,'ice density',True);finite(latent_J_kg,'latent heat',True)
    energy=net.diagnostics()['terminal_storage_to_ice_factor']*net.capacity_J_K*debris_temperature_excursion_K
    return dict(phase_energy_bound_J=energy,equivalent_lowering_bound_m=energy/(rho_i*latent_J_kg*ice_area_m2),
                scope='fixed network, supplied terminal temperature excursion; not a recovery guarantee')
