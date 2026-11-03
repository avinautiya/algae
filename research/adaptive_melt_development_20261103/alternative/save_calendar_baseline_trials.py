"""Preserve all predeclared strong-calendar baseline trials, DEVELOPMENTonly."""
from pathlib import Path
import json,hashlib
import numpy as np
from run_development import load_site,prediction_ratios,R3
P=Path(__file__).parent;out=P/'results_v1/CALENDAR_BASELINE_ALL_TRIALS.json'
if out.exists():raise ValueError('Do not overwrite baseline evidence')
r=json.loads((P/'results_v1/DEVELOPMENT_RESULTS.json').read_text());rows=[]
for tower in ['Tower1','Tower2']:
 f=load_site(tower);truth=np.array(r['truth'][tower]['development']);y0=r['truth'][tower]['train']
 for gamma in np.linspace(-.04,.04,801):
  p=y0*prediction_ratios(f,'calendar_exp_air',gamma)
  rows.append({'tower':tower,'gamma_per_day':float(gamma),'predictions_mmwe':p.tolist(),'absolute_development_errors_mmwe':abs(p-truth).tolist()})
with out.open('x')as file:json.dump({'scope':'All1602strongbaselineparametertrials retained; nofreshvalidation','source_anomaly_index_sha256':hashlib.sha256((R3/'moisture/anomaly_index.py').read_bytes()).hexdigest(),'trials':rows},file,indent=2);file.write('\n')
print('All1602calendar baseline trials preserved.')
