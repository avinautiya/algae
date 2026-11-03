"""Independent algebra and saved-result arithmetic; no campaign or fitted model."""
import json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def require_fixed_geometry(rows):
 """A convergence gate, not a generator of corrected masks or a new solver."""
 keys=['divider_width_m','total_pool_width_m','lane_width_m','pool_length_m','reach_length_m','tail_length_m']
 if len(rows)<2:raise ValueError('at least two meshes required')
 for k in keys:
  values=[r.get(k) for r in rows]
  if any(type(v) not in (int,float) or not math.isfinite(v) or v<=0 for v in values):raise ValueError('missing physical geometry '+k)
  if any(abs(v-values[0])>1e-12*max(1,abs(values[0])) for v in values[1:]):raise ValueError('not fixed-domain refinement: '+k)
 return True
RHOI,L,RHOW,CP=917.,334000.,1000.,4217.
Q0,T,H,GAP,SECONDS=.07,.5,200.,.25,21*86400
K=(2/3)**1.5*math.sqrt(9.81)*.5
HC=(Q0/K)**(2/3)
def melt(q,a,h=H):
 f=RHOW*CP*q
 return f*T*(-math.expm1(-h*a/f))/(RHOI*L*a)
def cap(q,h=H):return h*T/(RHOI*L)-melt(q/2,1.,h)
def split(q,g):
 lo=0.;hi=g+(q/K)**(2/3)
 for _ in range(80):
  mid=(lo+hi)/2
  if K*max(mid-g,0)**1.5+K*mid**1.5<q:lo=mid
  else:hi=mid
 eta=(lo+hi)/2;return K*max(eta-g,0)**1.5,K*eta**1.5
def small_area_repair(q=.14,area=.001):
 def rate(g):
  qt,qb=split(q,g);return melt(qb,1.)-melt(qt,area)
 g=GAP;t=0.;best=0.
 while t<SECONDS:
  dt=min(900.,SECONDS-t);g+=dt*rate(g+.5*dt*rate(g));t+=dt;best=max(best,GAP-g)
 return dict(final_gap_m=g,repair_m=best,reaches_critical_gap=g<=HC)
def root(fun,target,lo,hi):
 for _ in range(100):
  x=(lo+hi)/2
  if fun(x)<target:lo=x
  else:hi=x
 return (lo+hi)/2
def run():
 saved=json.loads((ROOT/'source/viability/viability_check.json').read_text());rows=[]
 for m in [1,2,4,8]:
  bound=cap(Q0*m)*SECONDS
  refs=[r for r in saved['cases'] if r['id'].endswith('_Q'+str(m))]
  assert all(abs(bound/r['bound_21d_m']-1)<1e-10 for r in refs)
  rows.append(dict(flow_multiple=m,bound_21days_m=bound,required_m=GAP-HC,viable=bound>=GAP-HC))
 assert not any(r['viable'] for r in rows)
 phi=H*T/(RHOI*L);n=H/(RHOW*CP*Q0/2)
 required_n=root(lambda x:1+math.expm1(-x)/x,(GAP-HC)/(phi*SECONDS),1e-8,10.)
 required_h=root(lambda h:cap(Q0,h)*SECONDS,GAP-HC,H,100000.)
 data=json.loads((ROOT/'source/verification/refinement.json').read_text())['controls']
 flat=[r for group in data.values() for r in group]
 maxima={k:max(r['max_over_steps'][k] for r in flat) for k in ['mass_kg','sensible_J','thermal_latent_J']}
 geom=[dict(dx_m=dx,lane_cells=round(.25/dx),divider_width_m=dx,total_pool_width_m=(2*round(.25/dx)+1)*dx) for dx in [.25,.125,.0625]]
 repair=small_area_repair()
 result=dict(audited_commit='455bf85856bcf11771f0b37bd63cc58e89f05f32',no_go_reproduced=True,cases=rows,nonzero_small_area_repair_counterexample=repair,report_all_area_sweep_repairs_zero_claim_is_false=True,declared_ntu=n,required_ntu_if_film_limit_fixed=required_n,ntu_ratio_if_film_limit_fixed=required_n/n,
  h_W_m2_K_removing_bound_if_h_itself_varied=required_h,h_ratio_if_h_itself_varied=required_h/H,ntu_when_h_itself_varied=required_h/(RHOW*CP*Q0/2),
  distinction='Both are necessary feasibility bounds only, not sufficient recovery or recommendations. NTU0.22 depends on keeping film-limit h*T fixed.',all_controls_per_step_maxima=maxima,refinement_geometry=geom,
  spatial_claim='These saved cases change physical pool/divider width. Differences are not pure fixed-domain mesh error. Toe partition remains unresolved.',can_numeric_or_wording_errors_reverse_declared_no_go=False)
 (ROOT/'CHECK_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
 return result
if __name__=='__main__':run()
