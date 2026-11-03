"""Conservative effective radiation-response reservoir; NOT complete glacier SEB.

Q=(1-a)SWIN+eps*(LWIN-sigma*Tburied^4), with fixed source a=.2,eps=.94.
E' = Q-max(E,0)/tau. Negative E is an EFFECTIVE deficit, not observed cold ice.
Positive outflow can be partitioned into calibrated latent heat and unobserved
loss; independent atmospheric H/LE/snow/ice-state terms are absent.
"""
import math
SIGMA=5.670374419e-8
LATENT=334000.
TIMESCALE_DAYS=(0.,.25,1.,3.)
ALBEDO=.2
EMISSIVITY=.94

def radiation_proxy(sw,lw,buried_C):
    if not all(math.isfinite(v) for v in [sw,lw,buried_C]) or sw<0 or lw<0 or buried_C<=-273.15:
        raise ValueError('Invalid radiation/temperature input')
    return (1-ALBEDO)*sw+EMISSIVITY*(lw-SIGMA*(buried_C+273.15)**4)

def step(energy,Q,seconds,tau):
    if not all(math.isfinite(v) for v in [energy,Q,seconds,tau]) or seconds<=0 or tau<0:
        raise ValueError('Finite state/positive step/nonnegative timescale required')
    initial=energy;left=seconds;delivered=0.
    if tau==0:
        delivered=max(energy,0.);energy=min(energy,0.)
        if Q>0:
            warm=min(left,-energy/Q);energy+=Q*warm;left-=warm
            if left>0:delivered+=Q*left;energy=0.
        else:energy+=Q*left
    else:
        if energy<=0:
            if Q<=0:return energy+Q*left,0.
            warm=min(left,-energy/Q);energy+=Q*warm;left-=warm
        if left>0:
            if Q<0:
                hit=tau*math.log1p(energy/(-Q*tau))
                active=min(left,hit)
            else:active=left
            loss=-math.expm1(-active/tau)
            new=energy*(1-loss)+Q*tau*loss
            if Q<0 and active==hit:new=0.
            delivered+=Q*active-(new-energy)
            energy=new;left-=active
            if left>0:energy+=Q*left
    residual=Q*seconds-(energy-initial)-delivered
    if abs(residual)>2e-9*max(1,abs(Q*seconds),abs(initial),abs(energy)):
        raise ArithmeticError('Effective reservoir ledger failed')
    if delivered < -1e-8:raise ArithmeticError('Negative delivered energy')
    return energy,delivered
