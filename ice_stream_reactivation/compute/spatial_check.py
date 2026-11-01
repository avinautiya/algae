from dataclasses import replace
from pathlib import Path
import json
from model import Case,evolve,partition
out=Path(__file__).parent/'records'
base=Case(discharge=.07,initial_debris_depth=.2,allow_supply=False)
cases=[]
for fraction in (.2,.4,.6,.8,1.):
 for footprint in (False,True):
  cases.append(replace(base,source_fractions=((1-fraction)/2,(1-fraction)/2,fraction),footprint_shield=footprint))
r=evolve(cases,days=45,dt=225)
(out/'spatial_coverage_check.json').write_text(json.dumps(r,indent=2))
for i,c in enumerate(cases):print(c.source_fractions[-1],c.footprint_shield,r['wet_days'][i],r['toe_melt_m3'][i],flush=True)
# Dense coarse deposit: no ongoing supply, no unbounded aggradation.
dense=replace(base,source_fractions=(.1,.1,.8))
rows=[]
for day in (1.,7.,14.,21.):
 for amp in (2.,4.,8.):
  p=dict(start_day=day,multiplier=amp,duration_day=.25)
  a=evolve([dense],days=45,dt=225,pulse=p,keep_history=False)
  K=(2/3)**1.5*dense.gravity**.5*dense.width
  q,_,_=partition(dense.discharge,a['final_toe_z'][0],a['final_bypass_z'][0],K,K)
  rows.append(dict(pulse=p,wet_days=a['wet_days'][0],toe_melt_m3=a['toe_melt_m3'][0],final_normal_toe_flow=float(q),final_stock=a['final_stock_m3'][0]))
(out/'spatial_pulse_check.json').write_text(json.dumps(rows,indent=2))
print('dense pulse restored',sum(a['final_normal_toe_flow']>1e-8 for a in rows),'of',len(rows),flush=True)
