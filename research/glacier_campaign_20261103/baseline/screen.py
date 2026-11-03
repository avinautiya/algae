"""Declared ideal two-route viability screen, not field or full-SWE prediction."""
import json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent
RHOI,L,RHOW,CP,H,T=917.,334000.,1000.,4217.,200.,.5
K=(2/3)**1.5*math.sqrt(9.81)*.5
Q0=.07; GAP=.25; HC=(Q0/K)**(2/3); HORIZON_D=21.
def rates(q,a):
 if q==0:return 0.
 f=RHOW*CP*q
 return f*T*(-math.expm1(-H*a/f))/(RHOI*L*a)
def partition(q,g):
 lo=0.;hi=g+(q/K)**(2/3)
 for _ in range(80):
  eta=(lo+hi)/2
  if K*max(eta-g,0)**1.5+K*eta**1.5<q:lo=eta
  else:hi=eta
 eta=(lo+hi)/2
 return K*max(eta-g,0)**1.5,K*eta**1.5
def run():
 rows=[]
 for at in [1.,.5,.25]:
  for mult in [1.,2.,4.,8.]:
   q=Q0*mult;qt,qb=partition(q,GAP)
   # For every g>=HC>0, equal controls imply qb>=q/2.
   # e(q/A) increases; toe <= asymptotic film heat limit.
   cap=H*T/(RHOI*L)-rates(q/2,1.)
   needed=GAP-HC
   rows.append(dict(id=f'A{at:g}_Q{mult:g}',toe_area_m2=at,bypass_area_m2=1.,pulse_Q_m3_s=q,initial_toe_Q_m3_s=qt,initial_repair_rate_m_day=(rates(qt,at)-rates(qb,1.))*86400,rigorous_repair_upper_bound_m_day=cap*86400,required_repair_m=needed,minimum_repair_days=needed/(cap*86400),maximum_possible_21day_repair_m=cap*86400*HORIZON_D,can_repair_in_21days_under_ideal_assumptions=cap*86400*HORIZON_D>needed,equal_area_no_reset_proof_applies=at==1.))
 return dict(scope='equal-K two-route controls, common fixed h/T, independent routes, no other heat or geometry sources; not universal environmental or SWE bound',horizon_days=HORIZON_D,source_parameters=dict(Q0_m3_s=Q0,K=K,initial_gap_m=GAP,critical_baseline_head_m=HC,h_W_m2_K=H,water_superheat_K=T),cases=rows,development_cases_already_inspected=True,production_gate='NO_GO for claiming durable rescue in these ideal cases; test other physics only if independently justified')
if __name__=='__main__':
 r=run();(ROOT/'RESULTS.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(cases=len(r['cases']),viable=sum(x['can_repair_in_21days_under_ideal_assumptions'] for x in r['cases']),scope=r['scope'])))
