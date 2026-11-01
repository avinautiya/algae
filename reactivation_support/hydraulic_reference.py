"""ODE-free reference for flat, wide, prismatic channels with a critical exit.

No entrance losses, contractions, jumps or arbitrary bed steps are modeled.
Useful for verification geometries satisfying those assumptions, not as truth
for the existing abrupt-step 2-D layout. Derivation in DIAGNOSTIC_FINDINGS.md.
"""
from dataclasses import dataclass
import math
from scipy.optimize import brentq


def _length_function(u):
    if abs(u)<1e-3:
        return sum(((13/3)**(k-1)-(4/3)**(k-1))/math.factorial(k)*u**k for k in range(2,10))
    return (3/13)*math.expm1((13/3)*u)-(3/4)*math.expm1((4/3)*u)


def upstream_head(unit_q,length,manning,g=9.81):
    for v in (unit_q,length,manning,g):
        if isinstance(v,bool) or not math.isfinite(v):raise ValueError('finite numeric input required')
    if min(unit_q,length,manning)<0 or g<=0:raise ValueError('invalid hydraulic parameters')
    if unit_q==0:return 0.
    hc=(unit_q*unit_q/g)**(1/3)
    if manning==0 or length==0:return 1.5*hc
    target=manning*manning*unit_q*unit_q*length/hc**(13/3)
    upper=1.
    while _length_function(upper)<target:upper*=2
    u=brentq(lambda u:_length_function(u)-target,0.,upper,xtol=1e-14,rtol=1e-14)
    h=hc*math.exp(u)
    return h+unit_q*unit_q/(2*g*h*h)


def unit_discharge(head,length,manning,g=9.81):
    if not math.isfinite(head):raise ValueError('finite head required')
    # Validate even when the lane is dry.
    upstream_head(0.,length,manning,g)
    if head<=0:return 0.
    ideal=math.sqrt(g)*(2*head/3)**1.5
    if manning==0 or length==0:return ideal
    x=brentq(lambda x:upstream_head(ideal*x,length,manning,g)-head,0.,1.,xtol=1e-14,rtol=1e-13)
    return ideal*x


@dataclass(frozen=True)
class Lane:
    z:float
    width:float
    length:float=2.
    manning:float=.02
    def __post_init__(self):
        if not all(math.isfinite(x) for x in (self.z,self.width,self.length,self.manning)) or self.width<=0 or min(self.length,self.manning)<0:
            raise ValueError('invalid lane')


def partition(Q,lanes,g=9.81):
    if not math.isfinite(Q) or Q<=0 or not lanes:raise ValueError('positive total discharge and lanes required')
    lo=min(l.z for l in lanes)
    hi=max(l.z+upstream_head(Q/l.width,l.length,l.manning,g) for l in lanes)
    def flows(eta):return [l.width*unit_discharge(eta-l.z,l.length,l.manning,g) for l in lanes]
    eta=brentq(lambda eta:sum(flows(eta))-Q,lo,hi,xtol=1e-13,rtol=1e-13)
    qs=flows(eta)
    return dict(stage=eta,discharges=qs,fractions=[q/Q for q in qs],mass_residual=sum(qs)-Q,
                model='wide horizontal prismatic lanes; critical exit; zero entrance loss')
