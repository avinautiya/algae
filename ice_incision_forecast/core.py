"""Exact quasi-steady thermal map and sharp regional heat-allocation bounds.

Assumptions: temperate ice, bed-only sensible-to-latent exchange, all meltwater
joins liquid at Tm, no other sources/storage, prescribed inlet flow, and positive
effective conductance. These are conditional bounds, not full CFD predictions.
"""
from bisect import bisect_right
from dataclasses import dataclass
from math import exp, expm1, isfinite, log


def finite(x,name,minimum=None,strict=False):
    if not isfinite(x) or (minimum is not None and (x<minimum or (strict and x==minimum))):
        raise ValueError(f"Invalid {name}")
    return x


@dataclass(frozen=True)
class ThermalMap:
    inlet_Q_m3_s: float
    inlet_excess_K: float
    rho_water: float=1000.0
    rho_ice: float=917.0
    cp: float=4217.0
    latent_heat: float=334000.0

    def __post_init__(self):
        for name,value in vars(self).items():
            finite(value,name,0.0,name!="inlet_excess_K")

    @property
    def capacity_W_K(self):
        return self.rho_water*self.cp*self.inlet_Q_m3_s

    @property
    def inlet_power_W(self):
        return self.capacity_W_K*self.inlet_excess_K

    @property
    def ste(self):
        return self.cp*self.inlet_excess_K/self.latent_heat

    def attenuation(self,conductance_W_K):
        finite(conductance_W_K,"conductance",0.0)
        g=conductance_W_K/self.capacity_W_K
        if g==0 or self.ste==0:
            return g
        lo,hi=g/(1+self.ste),g
        for _ in range(64):
            s=(lo+hi)/2
            value=s+self.ste*(s+expm1(-s))
            if value>g:
                hi=s
            else:
                lo=s
        return (lo+hi)/2

    def state(self,conductance_W_K):
        s=self.attenuation(conductance_W_K)
        power=self.inlet_power_W*exp(-s)
        removed=self.inlet_power_W*(-expm1(-s))
        q=self.inlet_Q_m3_s+removed/(self.rho_water*self.latent_heat)
        return dict(s=s,power_W=power,removed_W=removed,Q_m3_s=q,
                    temperature_excess_K=power/(self.rho_water*self.cp*q))

    def region_power(self,prefix_W_K,region_W_K):
        finite(prefix_W_K,"prefix conductance",0.0)
        finite(region_W_K,"region conductance",0.0)
        a=self.attenuation(prefix_W_K)
        b=self.attenuation(prefix_W_K+region_W_K)
        return self.inlet_power_W*exp(-a)*(-expm1(-(b-a)))

    def total_conductance_from_temperature(self,outlet_excess_K):
        """Exact historical/ideal calibration only, not a future prediction input."""
        finite(outlet_excess_K,"outlet temperature",0.0,True)
        if self.inlet_excess_K==0 or outlet_excess_K>self.inlet_excess_K:
            raise ValueError("Require 0 < outlet excess <= positive inlet excess")
        ratio=outlet_excess_K/self.inlet_excess_K
        p_ratio=ratio*(1+self.ste)/(1+ratio*self.ste)
        s=-log(p_ratio)
        return self.capacity_W_K*(s+self.ste*(s+expm1(-s)))


@dataclass(frozen=True)
class ConductancePrior:
    lower_W_K: tuple
    upper_W_K: tuple
    total_lower_W_K: float
    total_upper_W_K: float
    exact_prefix_W_K: tuple=()

    def __post_init__(self):
        object.__setattr__(self,"lower_W_K",tuple(self.lower_W_K))
        object.__setattr__(self,"upper_W_K",tuple(self.upper_W_K))
        object.__setattr__(self,"exact_prefix_W_K",tuple(tuple(x) for x in self.exact_prefix_W_K))
        if not self.lower_W_K or len(self.lower_W_K)!=len(self.upper_W_K):
            raise ValueError("Matched nonempty cell bounds required")
        for lo,hi in zip(self.lower_W_K,self.upper_W_K):
            finite(lo,"lower conductance",0.0); finite(hi,"upper conductance",lo)
        finite(self.total_lower_W_K,"total lower",0.0)
        finite(self.total_upper_W_K,"total upper",self.total_lower_W_K)
        indices=[x[0] for x in self.exact_prefix_W_K]
        if indices!=sorted(set(indices)) or any(type(k) is not int or not 0<k<len(self.lower_W_K) for k in indices):
            raise ValueError("Exact prefix indices must be unique increasing interior cell boundaries")
        for _,value in self.exact_prefix_W_K:
            finite(value,"exact prefix",0.0)
        self.blocks()  # validates the full feasible set

    def blocks(self):
        """Independent blocks defined by exact prefix anchors; last total is a band."""
        anchors=[(0,0.0)]+list(self.exact_prefix_W_K)
        result=[]
        for j,(left,offset) in enumerate(anchors):
            if j+1<len(anchors):
                right,value=anchors[j+1]; blo=bhi=value-offset
            else:
                right=len(self.lower_W_K)
                blo=self.total_lower_W_K-offset; bhi=self.total_upper_W_K-offset
            lower=max(sum(self.lower_W_K[left:right]),blo)
            upper=min(sum(self.upper_W_K[left:right]),bhi)
            if lower>upper+1e-10*max(1,abs(lower),abs(upper)):
                raise ValueError("Conductance bounds/prefix constraints are infeasible")
            if lower>upper:
                lower=upper
            result.append(dict(left=left,right=right,offset=offset,lower=lower,upper=upper))
        return result

    def feasible(self,profile,tolerance=1e-8):
        if len(profile)!=len(self.lower_W_K) or any(not isfinite(x) for x in profile):
            return False
        if any(x<lo-tolerance or x>hi+tolerance for x,lo,hi in zip(profile,self.lower_W_K,self.upper_W_K)):
            return False
        total=sum(profile)
        if total<self.total_lower_W_K-tolerance or total>self.total_upper_W_K+tolerance:
            return False
        return all(abs(sum(profile[:k])-value)<=tolerance for k,value in self.exact_prefix_W_K)


def fill(prior,profile,left,right,target):
    lo=sum(prior.lower_W_K[left:right]); hi=sum(prior.upper_W_K[left:right])
    tol=1e-9*max(1,abs(target),abs(hi))
    if target<lo-tol or target>hi+tol:
        raise ValueError("Internal infeasible allocation")
    fraction=(target-lo)/(hi-lo) if hi>lo else 0.0
    fraction=max(0.0,min(1.0,fraction))
    for i in range(left,right):
        profile[i]=prior.lower_W_K[i]+fraction*(prior.upper_W_K[i]-prior.lower_W_K[i])


def default_profile(prior):
    out=[0.0]*len(prior.lower_W_K)
    for block in prior.blocks():
        fill(prior,out,block["left"],block["right"],(block["lower"]+block["upper"])/2)
    return out


def region_bounds(prior,start,stop):
    """Sharp regional bounds as conductance witnesses, valid for every ThermalMap.

Indices select contiguous cells [start,stop). Profiles attain the lower/upper
regional heat for ANY positive inlet flow/temperature and common nonnegative
time multiplier. Marginal bounds for different regions need different witnesses.
"""
    n=len(prior.lower_W_K)
    if type(start) is not int or type(stop) is not int or not 0<=start<stop<=n:
        raise ValueError("Require a nonempty contiguous region")
    blocks=prior.blocks(); lo=prior.lower_W_K; hi=prior.upper_W_K
    anchors={0:0.0,**{k:v for k,v in prior.exact_prefix_W_K}}
    def owner(k):
        if k in anchors:
            return None
        return next(i for i,b in enumerate(blocks) if b["left"]<k<=b["right"])
    ia,ib=owner(start),owner(stop)
    def endpoint(profile,k,index,want_low):
        if index is None:
            return anchors[k]
        b=blocks[index]; l,r=b["left"],b["right"]
        amin=max(sum(lo[l:k]),b["lower"]-sum(hi[k:r]))
        amax=min(sum(hi[l:k]),b["upper"]-sum(lo[k:r]))
        a=amin if want_low else amax
        tail=max(sum(lo[k:r]),b["lower"]-a)
        fill(prior,profile,l,k,a); fill(prior,profile,k,r,tail)
        return b["offset"]+a
    witnesses={}
    for label,maximize in (("minimum",False),("maximum",True)):
        profile=default_profile(prior)
        if ia is not None and ia==ib:
            block=blocks[ia]; l,r=block["left"],block["right"]
            al,au=sum(lo[l:start]),sum(hi[l:start])
            bl,bu=sum(lo[start:stop]),sum(hi[start:stop])
            cl,cu=sum(lo[stop:r]),sum(hi[stop:r])
            gl,gu=block["lower"],block["upper"]
            if maximize:
                a=max(al,gl-bu-cu)
                b=min(bu,gu-a-cl)
            else:
                a=min(au,gu-bl-cl)
                b=max(bl,gl-a-cu)
            c=max(cl,gl-a-b)
            fill(prior,profile,l,start,a); fill(prior,profile,start,stop,b); fill(prior,profile,stop,r,c)
        else:
            endpoint(profile,start,ia,want_low=maximize)
            endpoint(profile,stop,ib,want_low=not maximize)
        if not prior.feasible(profile,tolerance=1e-7*max(1,sum(prior.upper_W_K))):
            raise ArithmeticError("Constructed witness failed feasibility")
        witnesses[label]=dict(profile_W_K=profile,prefix_W_K=sum(profile[:start]),
            region_W_K=sum(profile[start:stop]),total_W_K=sum(profile))
    return dict(start=start,stop=stop,**witnesses,
                scope="Sharp conditional regional bounds; different regions have different attaining profiles")
