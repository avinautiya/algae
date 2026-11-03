"""Idealized linear cold-reference identities; no field fitting or melt scoring.

These are standard diffusion/admittance calculations, not originality claims.
"""
import cmath,math

def ice_impedance(s,conductivity,capacity,depth=None):
 if not all(type(v) in (int,float) and math.isfinite(v) and v>0 for v in [conductivity,capacity]):raise ValueError('positive ice material properties required')
 s=complex(s)
 if s.real<0 or not math.isfinite(s.real) or not math.isfinite(s.imag) or s==0:raise ValueError('nonzero physical Laplace/Fourier parameter required')
 root=cmath.sqrt(s);impedance=1/(math.sqrt(conductivity*capacity)*root)
 if depth is not None:
  if not math.isfinite(depth) or depth<=0:raise ValueError('positive ice depth required')
  impedance*=cmath.tanh(depth*cmath.sqrt(s*capacity/conductivity))
 return impedance

def length_response(s,k,contact_R,ice_k,ice_C,ice_depth=None):
 return k*(contact_R+ice_impedance(s,ice_k,ice_C,ice_depth))

def recover_two_real_laplace(s1,s2,ell1,ell2,alpha,ice_k,ice_C,ice_depth=None):
 if min(s1,s2,alpha)<=0 or s1==s2:raise ValueError('distinct positive real frequencies and positive diffusivity required')
 z1=ice_impedance(s1,ice_k,ice_C,ice_depth).real;z2=ice_impedance(s2,ice_k,ice_C,ice_depth).real
 if z1==z2:raise ValueError('reference impedances do not separate')
 k=(ell1-ell2)/(z1-z2);R=ell1/k-z1
 if k<=0 or R< -1e-12:raise ValueError('data inconsistent with positive conductivity/passive contact')
 return dict(conductivity=k,contact_R=R,capacity=k/alpha,
  gain_from_two_ell_errors_to_k=1/abs(z1-z2),
  scope='noise-free structural algebra only; assumes no melting and known initial/reference response',field_identifiability_established=False,mathematical_novelty_claim=False)

def warm_gauge(k,C,R,scale):
 return dict(alpha=k/C,contact_length=k*R,scaled_alpha=(scale*k)/(scale*C),scaled_contact_length=(scale*k)*(R/scale),melt_flux_scale=scale)
