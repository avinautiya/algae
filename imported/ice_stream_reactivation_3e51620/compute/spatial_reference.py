from pathlib import Path
import json
from model import Case,evolve
from independent_reference import reference
base=Case(discharge=.07,initial_debris_depth=.2,allow_supply=False,source_fractions=(.1,.1,.8),footprint_shield=True)
rows=[]
for pulse in (None,dict(start_day=1.,multiplier=2.,duration_day=.25),dict(start_day=14.,multiplier=8.,duration_day=.25)):
 ref=reference(base,days=45,pulse=pulse)
 fixed=evolve([base],days=45,dt=112.5,pulse=pulse,keep_history=False)
 rows.append(dict(pulse=pulse,reference=ref,fixed=fixed))
 print('comparison',pulse,ref['toe_melt_m3'],fixed['toe_melt_m3'][0],ref['events_days'],fixed['first_dry_day'],flush=True)
(Path(__file__).parent/'records/spatial_independent_reference.json').write_text(json.dumps(rows,indent=2))
