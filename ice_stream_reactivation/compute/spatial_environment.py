from dataclasses import replace
from pathlib import Path
import json
from model import Case,evolve
base=Case(discharge=.07,initial_debris_depth=.2,allow_supply=False,source_fractions=(.1,.1,.8),footprint_shield=True)
cases=[]
for temp in (.2,.5,1.):
 for h in (100.,200.,500.):
  for extra in (0.,25.,50.,100.,200.):
   cases.append(replace(base,water_temperature=temp,h_water=h,extra_toe_heat_W_m2=extra))
r=evolve(cases,days=45,dt=225,keep_history=False)
(Path(__file__).parent/'records/spatial_environment_sensitivity.json').write_text(json.dumps(r,indent=2))
print('synthetic sensitivity cases',len(cases),'contact loss',sum(x<45 for x in r['wet_days']),flush=True)
for i,c in enumerate(cases):
 if c.h_water==200 and c.water_temperature==.5:print(c.extra_toe_heat_W_m2,r['wet_days'][i],flush=True)
