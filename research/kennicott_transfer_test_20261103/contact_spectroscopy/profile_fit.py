"""Geometry-explicit extension; original 2020 implementation/results stay pinned."""
import math
from .spectroscopy import transfer

def _geometry(L,x):
 if not math.isfinite(L) or not math.isfinite(x) or L<=0 or x<=0 or x>L+1e-12:raise ValueError('Deep probe lies outside modeled debris geometry')
 return L,min(x,L)
def components(alpha,b,observations,L,x):
 L,x=_geometry(L,x)
 return [abs(transfer(alpha,b,z['omega_rad_s'],L,x)-complex(z['ratio_real'],z['ratio_imag']))/z['ratio_error_scale'] for z in observations]
def fit_geometry(obs,*,L,x,null=False):
 L,x=_geometry(L,x);amin,amax=math.log(1e-8),math.log(5e-6);lo,hi=amin,amax;blo,bhi=0.,.3;best=None
 for level in range(4):
  na=61 if level==0 else 21;nb=1 if null else na;da=(hi-lo)/(na-1);db=0 if null else (bhi-blo)/(nb-1)
  for i in range(na):
   a=lo+i*da
   for j in range(nb):
    b=0. if null else blo+j*db;score=sum(y*y for y in components(math.exp(a),b,obs,L,x));candidate=(score,a,b)
    if best is None or candidate<best:best=candidate
  lo=max(amin,best[1]-da);hi=min(amax,best[1]+da)
  if not null:blo=max(0.,best[2]-db);bhi=min(.3,best[2]+db)
 score,loga,b=best;a=math.exp(loga)
 return {'alpha_m2_s':a,'b_kR_m':b,'weighted_objective':score,'per_ratio_error_scale_residuals':components(a,b,obs,L,x),'alpha_at_bound':min(abs(loga-amin),abs(loga-amax))<1e-6,'b_at_bound':null or min(b,.3-b)<1e-7}
def rank_geometry(alpha,b,obs,*,L,x,h=1e-4):
 L,x=_geometry(L,x)
 def vec(a,bb):
  vals=[]
  for row in obs:
   t=transfer(a,bb,row['omega_rad_s'],L,x)/row['ratio_error_scale'];vals.append(t.real)
   if row['omega_rad_s']:vals.append(t.imag)
  return vals
 v=vec(alpha*math.exp(h),b);w=vec(alpha*math.exp(-h),b);j1=[(a-z)/(2*h) for a,z in zip(v,w)]
 bb=max(0,b-h*L);cc=b+h*L;v=vec(alpha,cc);w=vec(alpha,bb);j2=[(a-z)/((cc-bb)/L) for a,z in zip(v,w)]
 a=sum(t*t for t in j1);d=sum(t*t for t in j2);c=sum(s*t for s,t in zip(j1,j2));disc=math.sqrt((a-d)**2+4*c*c);s1=math.sqrt((a+d+disc)/2);s2=math.sqrt(max(0,(a+d-disc)/2));return {'coordinates':['logalpha','b/L'],'singular_values':[s1,s2],'condition_number':s1/s2 if s2 else None,'scope':'local fixed-geometry response, not global physical identifiability'}
