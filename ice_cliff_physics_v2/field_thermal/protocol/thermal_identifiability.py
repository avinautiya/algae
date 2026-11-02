"""Independent harmonic diffusion reference and flux-scale identifiability check.

This is established linear heat-equation mathematics, not a new theorem.
"""
import math,cmath

def _check(L,alpha,period):
    if any(isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<=0 for x in [L,alpha,period]):raise ValueError('positive finite L,alpha,period')

def _inverse_sinh(x):
    return 2*cmath.exp(-x)/(1-cmath.exp(-2*x))

def harmonic_temperature(z,t,*,length_m,alpha_m2_s,mean_top_C=5.,amplitude_top_C=2.,period_s=86400.):
    L=length_m;_check(L,alpha_m2_s,period_s)
    if not 0<=z<=L:raise ValueError('depth outside domain')
    gamma=cmath.sqrt(1j*2*math.pi/(period_s*alpha_m2_s))
    ratio=cmath.exp(-gamma*z)*(1-cmath.exp(-2*gamma*(L-z)))/(1-cmath.exp(-2*gamma*L))
    return mean_top_C*(1-z/L)+(amplitude_top_C*ratio*cmath.exp(1j*2*math.pi*t/period_s)).real

def harmonic_bottom_flux(t,*,length_m,alpha_m2_s,conductivity_W_m_K,mean_top_C=5.,amplitude_top_C=2.,period_s=86400.):
    _check(length_m,alpha_m2_s,period_s)
    if conductivity_W_m_K<=0:raise ValueError('positive conductivity')
    gamma=cmath.sqrt(1j*2*math.pi/(period_s*alpha_m2_s))
    return conductivity_W_m_K*(mean_top_C/length_m+(amplitude_top_C*gamma*_inverse_sinh(gamma*length_m)*cmath.exp(1j*2*math.pi*t/period_s)).real)

def harmonic_weighted_moment(t,*,length_m,alpha_m2_s,mean_top_C=5.,amplitude_top_C=2.,period_s=86400.):
    L=length_m;_check(L,alpha_m2_s,period_s)
    gamma=cmath.sqrt(1j*2*math.pi/(period_s*alpha_m2_s))
    factor=1/(L*gamma*gamma)-_inverse_sinh(gamma*L)/gamma
    return mean_top_C*L/6+(amplitude_top_C*factor*cmath.exp(1j*2*math.pi*t/period_s)).real
