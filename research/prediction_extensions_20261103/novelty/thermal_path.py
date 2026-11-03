"""Classical serial heat exchange; synthetic counterexample, not a new law."""
from dataclasses import dataclass,asdict
import math,json

CW=1000*4217.
LAT=917*334000.
@dataclass(frozen=True)
class Path:
 flow_m3_s:float
 upstream_conductance_W_K:float
 target_conductance_W_K:float
 target_area_m2:float
 inlet_superheat_K:float
def evaluate(p):
 for k,v in asdict(p).items():
  if type(v) not in (float,int) or not math.isfinite(v) or v<0:raise ValueError(k)
 if p.target_area_m2<=0:raise ValueError('positive target area required')
 f=CW*p.flow_m3_s;power=f*p.inlet_superheat_K
 if f==0:return dict(inlet_W=0.,upstream_melt_W=0.,target_melt_W=0.,outlet_W=0.,target_inlet_K=None,incision_m_s=0.,residual_W=0.)
 u=math.exp(-p.upstream_conductance_W_K/f)
 t=math.exp(-p.target_conductance_W_K/f)
 up=power*(-math.expm1(-p.upstream_conductance_W_K/f))
 target=power*u*(-math.expm1(-p.target_conductance_W_K/f))
 out=power*u*t
 return dict(inlet_W=power,upstream_melt_W=up,target_melt_W=target,outlet_W=out,target_inlet_K=p.inlet_superheat_K*u,incision_m_s=target/(LAT*p.target_area_m2),residual_W=power-up-target-out)
def example():
 # Hydraulics are an audited ideal-control pulse illustration, not measurements.
 qt=.006834176791590543;qb=.13316582320840942
 a=Path(qt,0.,100.,.5,.5);b=Path(qb,0.,200.,1.,.5)
 x=evaluate(a);y=evaluate(b)
 # Solve exact equal-incision upstream-bypass conductance; no target fitting.
 critical=CW*qb*math.log(y['incision_m_s']/x['incision_m_s'])
 cooled=Path(qb,2*critical,200.,1.,.5);z=evaluate(cooled)
 return dict(scope='synthetic prescribed-flow/path example; no field validation or novelty claim',toe=asdict(a),bypass_direct=asdict(b),bypass_longer=asdict(cooled),critical_upstream_bypass_conductance_W_K=critical,
  no_upstream_paths=dict(toe=x,bypass=y,gap_drift_m_day=(y['incision_m_s']-x['incision_m_s'])*86400),
  upstream_bypass_path=dict(toe=x,bypass=z,gap_drift_m_day=(z['incision_m_s']-x['incision_m_s'])*86400),
  interpretation='Same target hydraulics; route-specific prior heat loss reverses local incision ranking. Heat deposited upstream remains glacier melt and cannot be counted as lost energy.')
if __name__=='__main__':print(json.dumps(example(),indent=2,allow_nan=False))
