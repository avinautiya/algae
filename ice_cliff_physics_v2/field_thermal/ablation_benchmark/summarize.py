import json
from pathlib import Path
import numpy as np
p=Path(__file__).parent;r=json.loads((p/'results.json').read_text());rows=r['source_target_comparisons'];out={}
for geom in ['profile','stake']:
 out[geom]={}
 for fam in ['pde','quasi_steady','one_node','constant_daily']:
  primary=[a for a in rows if a['geometry']==geom and a['family']==fam and a['source']=='2020b' and a['start_shift_h']==a['end_shift_h']==0 and a['early_geometry']=='current']
  allprimary=[a for a in rows if a['geometry']==geom and a['family']==fam and a['source']=='2020b']
  allfold=[a for a in rows if a['geometry']==geom and a['family']==fam and a['start_shift_h']==a['end_shift_h']==0 and a['early_geometry']=='current']
  out[geom][fam]={'primary_default_clock_current_early':[{'target':a['target'],'predicted_kg_m2':a['predicted_full_kg_m2'],'observed_kg_m2':a['observed_kg_m2'],'relative_error':a['relative_error']} for a in primary], 'primary_all_scenario_relative_error_minmax':[min(a['relative_error'] for a in allprimary),max(a['relative_error'] for a in allprimary)],'primary_default_MAPE':float(np.mean([abs(a['relative_error']) for a in primary])),'all12_default_MAPE':float(np.mean([abs(a['relative_error']) for a in allfold])),'all12_default_MAE_kg_m2':float(np.mean([abs(a['predicted_full_kg_m2']-a['observed_kg_m2']) for a in allfold])),'covered_exceeds_full_all_primary_scenarios':sum(a.get('covered_exceeds_full',False) for a in allprimary)}
(p/'SUMMARY.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
