"""Classical eigenfunction reference for a zero-initial one-dimensional slab.

No numerical diffusion solver or lumped ODE is imported. This is an analytic
reference implementation, not a new theorem or an environmental validation.
"""
from __future__ import annotations
import math
import sys


def _number(value,name,positive=False,nonnegative=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError('finite '+name+' required')
    if positive and value<=0 or nonnegative and value<0:raise ValueError('invalid '+name+' sign')
    return float(value)


def water_exchange_conductance(water_capacity_W_K,debris_water_conductance_W_K):
    F=_number(water_capacity_W_K,'water capacity',nonnegative=True)
    G=_number(debris_water_conductance_W_K,'debris conductance',nonnegative=True)
    return F*(-math.expm1(-G/F)) if F>0 else 0.


def eigenvalues(B,terms):
    """lambda cos(lambda)+B sin(lambda)=0, monotone phase-equation roots."""
    _number(B,'Robin ratio',nonnegative=True)
    if isinstance(terms,bool) or not isinstance(terms,int) or not 1<=terms<=4096:raise ValueError('terms must be an integer from1 to4096')
    out=[]
    for n in range(terms):
        base=(n+.5)*math.pi
        if B==0:lam=base
        else:
            lo=base;hi=(n+1)*math.pi
            for _ in range(64):
                mid=(lo+hi)/2
                if mid==lo or mid==hi:break
                value=mid-base-math.atan2(B,mid)
                if value<0:lo=mid
                else:hi=mid
            lam=(lo+hi)/2
        out.append(lam)
    return tuple(out)


class RobinSlab:
    """Zero-initial slab with bottom0 C and K*T_xi(1)=P-H*T(1).

    C and K are total capacity[J/K] and through-slab conductance[W/K].
    For physical water exchange, H=F*(1-exp(-Gd/F)), P=Qa+H*Tin.
    Use from_linear_boundary for a generic affine boundary, such as a separately
    derived parallel-contact q_debris=a*Tin-b*Ts; its top heat is not automatically
    the entire stream heat. Coefficients/forcing are constant over the solution.
    """
    def __init__(self,capacity_J_K,ice_conductance_W_K,water_exchange_W_K,
                 inlet_temperature_C,atmospheric_power_W=0.,terms=256):
        self.C=_number(capacity_J_K,'capacity',positive=True)
        self.K=_number(ice_conductance_W_K,'ice conductance',positive=True)
        self.H=_number(water_exchange_W_K,'Robin loss',nonnegative=True)
        self.Tin=_number(inlet_temperature_C,'inlet temperature',nonnegative=True)
        self.Qa=_number(atmospheric_power_W,'atmospheric power')
        self.P=self.Qa+self.H*self.Tin
        if not math.isfinite(self.P) or self.P<0:raise ValueError('finite nonnegative top forcing required')
        self.tau=self.C/self.K;self.B=self.H/self.K
        if not math.isfinite(self.tau) or self.tau<=0 or not math.isfinite(self.B):raise ValueError('unrepresentable diffusion scale or Robin ratio')
        if not math.isfinite(self.K+self.H):raise ValueError('conductance sum overflow')
        self.Ts=self.P/(self.K+self.H)
        self.terms=terms;self.roots=eigenvalues(self.B,terms);self.source_interpretation='water plus atmospheric boundary'
        coeff=[];sines=[];cosines=[];means=[]
        for n,lam in enumerate(self.roots):
            radius=math.hypot(lam,self.B);sign=1 if n%2==0 else -1
            sine=sign*lam/radius;cosine=-sign*self.B/radius
            norm=.5-sine*cosine/(2*lam)
            projection=-cosine/lam+sine/(lam*lam)
            coeff.append(-self.Ts*projection/norm);sines.append(sine);cosines.append(cosine)
            one_minus_cos=sine*sine/(1+cosine) if cosine>0 else 1-cosine
            means.append(one_minus_cos/lam)
        self.coefficients=tuple(coeff);self.sines=tuple(sines);self.cosines=tuple(cosines);self.spatial_means=tuple(means)

    @classmethod
    def from_contact(cls,capacity_J_K,ice_conductance_W_K,water_capacity_W_K,
                     debris_water_conductance_W_K,inlet_temperature_C,
                     atmospheric_power_W=0.,terms=256):
        return cls(capacity_J_K,ice_conductance_W_K,
                   water_exchange_conductance(water_capacity_W_K,debris_water_conductance_W_K),
                   inlet_temperature_C,atmospheric_power_W,terms)

    @classmethod
    def from_linear_boundary(cls,capacity_J_K,ice_conductance_W_K,
                             source_power_W,loss_conductance_W_K,terms=256):
        obj=cls(capacity_J_K,ice_conductance_W_K,loss_conductance_W_K,0.,source_power_W,terms)
        obj.source_interpretation='generic affine boundary P-H*Ts; no physical source split asserted'
        return obj

    def _time(self,time_s):return _number(time_s,'time',nonnegative=True)

    def _xi(self,xi):
        _number(xi,'normal coordinate')
        if not 0<=xi<=1:raise ValueError('normal coordinate must lie in[0,1]')
        return xi

    def _sine(self,lam,sine,xi):
        return 0. if xi==0 else sine if xi==1 else math.sin(lam*xi)

    def temperature(self,xi,time_s):
        xi=self._xi(xi);t=self._time(time_s)
        if t==0:return 0.
        return self.Ts*xi+math.fsum(c*self._sine(lam,s,xi)*math.exp(-lam*lam*t/self.tau)
            for lam,c,s in zip(self.roots,self.coefficients,self.sines))

    def initial_series_residual(self,xi):
        """Un-special-cased truncated initial series, useful for convergence QA."""
        xi=self._xi(xi)
        return self.Ts*xi+math.fsum(c*self._sine(lam,s,xi) for lam,c,s in zip(self.roots,self.coefficients,self.sines))

    def spatial_derivative(self,xi,time_s):
        xi=self._xi(xi);t=self._time(time_s)
        if t==0:raise ValueError('forced boundary is not classically compatible at initial time')
        return self.Ts+math.fsum(c*lam*(co if xi==1 else math.cos(lam*xi))*math.exp(-lam*lam*t/self.tau)
            for lam,c,co in zip(self.roots,self.coefficients,self.cosines))

    def cell_average(self,low_xi,high_xi,time_s):
        lo=self._xi(low_xi);hi=self._xi(high_xi);t=self._time(time_s)
        if hi<=lo:raise ValueError('positive cell width required')
        if t==0:return 0.
        mid=(lo+hi)/2;width=hi-lo
        def mean_sine(lam):
            arg=lam*width/2
            return math.sin(lam*mid)*(math.sin(arg)/arg if arg else 1.)
        return self.Ts*mid+math.fsum(c*mean_sine(lam)*math.exp(-lam*lam*t/self.tau) for lam,c in zip(self.roots,self.coefficients))

    def response(self,time_s):
        t=self._time(time_s)
        if t==0:
            return dict(time_s=0.,ice_power_W=0.,ice_heat_J=0.,top_temperature_C=0.,top_temperature_integral_K_s=0.,
                        water_heat_J=0. if self.source_interpretation.startswith('water') else None,
                        atmospheric_heat_J=0. if self.source_interpretation.startswith('water') else None,
                        boundary_heat_J=0.,stored_energy_J=0.,energy_residual_J=0.,roundoff_cancellation_scale_J=0.,
                        terms=self.terms,source_interpretation=self.source_interpretation)
        decay=[math.exp(-lam*lam*t/self.tau) for lam in self.roots]
        # Exact sum identities remove algebraic tails from integrated series.
        a=(1+2/(1+self.B))/6
        surface_lag=(1/(1+self.B))/3
        ice_exp=math.fsum(c/lam*e for lam,c,e in zip(self.roots,self.coefficients,decay))
        top_exp=math.fsum(c*s/(lam*lam)*e for lam,c,s,e in zip(self.roots,self.coefficients,self.sines,decay))
        ice=self.K*(self.Ts*(t-self.tau*a)-self.tau*ice_exp)
        top_integral=self.Ts*(t-self.tau*surface_lag)-self.tau*top_exp
        stored=self.C*(self.Ts/2+math.fsum(c*mean*e for c,mean,e in zip(self.coefficients,self.spatial_means,decay)))
        power=self.K*(self.Ts+math.fsum(c*lam*e for lam,c,e in zip(self.roots,self.coefficients,decay)))
        top=self.Ts+math.fsum(c*s*e for c,s,e in zip(self.coefficients,self.sines,decay))
        boundary=self.P*t-self.H*top_integral
        water=self.H*(self.Tin*t-top_integral) if self.source_interpretation.startswith('water') else None
        air=self.Qa*t if water is not None else None
        residual=(air+water if water is not None else boundary)-ice-stored
        absolute_work=(abs(self.P*t)+abs(self.H*top_integral)+abs(ice)+abs(stored)+
            self.K*(abs(self.Ts*(t-self.tau*a))+self.tau*math.fsum(abs(c/lam*e) for lam,c,e in zip(self.roots,self.coefficients,decay))))
        roundoff_scale=64*sys.float_info.epsilon*absolute_work
        values=[ice,top_integral,stored,power,top,boundary,residual,roundoff_scale]
        if not all(math.isfinite(x) for x in values):raise ValueError('derived analytic inventory overflow')
        return dict(time_s=t,ice_power_W=power,ice_heat_J=ice,top_temperature_C=top,
            top_temperature_integral_K_s=top_integral,water_heat_J=water,atmospheric_heat_J=air,
            boundary_heat_J=boundary,stored_energy_J=stored,energy_residual_J=residual,
            roundoff_cancellation_scale_J=roundoff_scale,terms=self.terms,source_interpretation=self.source_interpretation)
