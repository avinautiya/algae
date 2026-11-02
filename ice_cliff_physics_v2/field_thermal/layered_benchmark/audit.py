import csv,json,hashlib
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent;F=P.parent
R=json.loads((P/'results/RESULTS.json').read_text());checks=[];maxerr=0
for r in R['profiles']:
 name=r['profile_id'];p=P/'results'/f'{name}.csv';assert hashlib.sha256(p.read_bytes()).hexdigest()==r['prediction_sha256']
 rows=list(csv.DictReader(p.open()));source=list(csv.DictReader((F/'source/Debris_Temperature_Profiles/Profile_Temperature_Records'/f'{name}.csv').open(encoding='utf-8-sig')));lookup={x['Datetime']:x for x in source}
 depths=sorted({float(x['depth_m']) for x in rows});n=0
 for row in rows:
  z=float(row['depth_m']);raw=lookup[row['utc']];key=next(k for k in raw if k.startswith('Depth:') and abs(float(k.split(':')[1].split()[0])/100-z)<1e-9);assert float(raw[key])==float(row['observed_C']);n+=1
 for m,metric in r['metrics'].items():
  maes=[];mses=[]
  for z in depths:
   es=np.array([float(x[m])-float(x['observed_C']) for x in rows if float(x['depth_m'])==z]);maes.append(float(np.mean(abs(es))));mses.append(float(np.mean(es**2)))
  err=max(abs(np.mean(maes)-metric['MAE_C']),abs(np.mean(mses)-metric['MSE_C2']));maxerr=max(maxerr,err);assert err<1e-10
 checks.append({'profile':name,'raw_truth_rows_matched':n,'metrics_pass':True})
(P/'results/INDEPENDENT_AUDIT.json').write_text(json.dumps({'pass':True,'checks':checks,'max_metric_disagreement':maxerr},indent=2)+'\n');print('audit passed',sum(x['raw_truth_rows_matched'] for x in checks),'rows',maxerr)
