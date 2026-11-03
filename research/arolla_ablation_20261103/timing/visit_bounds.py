"""Exact visit-time bounds for one-driver calibrated interval forecasts.

Established monotonicity, not a new mathematical priority claim. Input drivers
must be nonnegative on explicit fully covered time bins. No weather is filled.
"""
import math
def finite(x):
 try:return type(x) in (float,int) and math.isfinite(x)
 except OverflowError:return False
def integrate(segments,start,end):
 if not finite(start) or not finite(end) or not start<end:raise ValueError('ordered finite interval required')
 cursor=start;total=0.;last=None
 for a,b,v in segments:
  if not all(finite(x) for x in [a,b,v]) or a>=b or v<0 or (last is not None and a<last):raise ValueError('ordered nonnegative finite driver bins required')
  last=b
  lo=max(a,start);hi=min(b,end)
  if hi<=lo:continue
  if abs(lo-cursor)>1e-8:raise ValueError('unobserved driver support; no extrapolation')
  total+=(hi-lo)*v;cursor=hi
 if abs(cursor-end)>1e-8 or not math.isfinite(total):raise ValueError('incomplete or overflowing driver integral')
 return total
def ratio(segments,t0,t1,t2):
 den=integrate(segments,t0,t1)
 if den<=0:raise ValueError('training coefficient unidentified')
 return integrate(segments,t1,t2)/den
def exact_bounds(segments,first,middle,last):
 if any(not isinstance(x,(tuple,list)) or len(x)!=2 or not all(finite(v) for v in x) or x[0]>x[1] for x in [first,middle,last]):raise ValueError('finite ordered visit ranges required')
 if not first[1]<middle[0] or not middle[1]<last[0]:raise ValueError('visit uncertainty ranges overlap')
 integrate(segments,first[0],last[1])
 lo=ratio(segments,first[0],middle[1],last[0]);hi=ratio(segments,first[1],middle[0],last[1])
 return dict(minimum_ratio=lo,maximum_ratio=hi,shared_middle_visit=True,interpretation='exact deterministic bounds for the supplied nonnegative driver and time ranges; not a confidence interval')
def residual_uncertainty_bound(factor,delta,*,association):
 if not finite(factor) or factor<0 or not finite(delta) or delta<0:raise ValueError('nonnegative finite factor and bounded error required')
 if association=='each_shared_reading':return 2*(factor+1)*delta
 if association=='each_interval':return (factor+1)*delta
 raise ValueError('source uncertainty association is not known; specify a labelled sensitivity')
