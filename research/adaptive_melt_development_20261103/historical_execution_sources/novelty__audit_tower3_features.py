"""Feature-only spatial transfer audit. Never imports camera values."""
import csv,json,math,bisect,hashlib
from datetime import datetime
from pathlib import Path
R=Path('/private/tmp/glacier-round2-20261103/data/decoded')
def load(n):
 p=R/n
 with p.open() as f:r=list(csv.DictReader(f))
 return [(datetime.fromisoformat(x['SOURCE_LOCAL_LABEL']).timestamp(),float(x['TA_3'])) for x in r],hashlib.sha256(p.read_bytes()).hexdigest()
def interp(rows,t):
 ts=[x[0] for x in rows];i=bisect.bisect_left(ts,t)
 if i<len(ts) and ts[i]==t:return rows[i][1],0
 if i==0 or i==len(ts):return None,None
 a,b=rows[i-1],rows[i];w=(t-a[0])/(b[0]-a[0]);return (1-w)*a[1]+w*b[1],b[0]-a[0]
a,ha=load('trh_Tower1.csv');b,hb=load('trh_Tower3.csv')
pairs=[]
for t,y in b:
 x,g=interp(a,t)
 if x is not None and math.isfinite(x) and math.isfinite(y) and g<=18*3600:pairs.append((t,x,y,g))
split=(pairs[0][0]+pairs[-1][0])/2
tr=[r for r in pairs if r[0]<split];te=[r for r in pairs if r[0]>=split]
xm=sum(r[1] for r in tr)/len(tr);ym=sum(r[2] for r in tr)/len(tr)
beta=sum((r[1]-xm)*(r[2]-ym) for r in tr)/sum((r[1]-xm)**2 for r in tr);alpha=ym-beta*xm
errors=[alpha+beta*r[1]-r[2] for r in te]
out={'scope':'Feature-only audit, no camera values read; not proof of future climate stationarity or melt accuracy.', 'source_hashes':{'Tower1':ha,'Tower3':hb},'clock':'Native naive local labels; timestamp numbers used only as same-clock differences, not inferred UTC.', 'method':'Tower1 TA_3 chosen standing tower, affine fit on chronological first half, second half held out; T1 interpolation <=18h.', 'train_n':len(tr),'test_n':len(te),'first':datetime.fromtimestamp(pairs[0][0]).isoformat(),'split':datetime.fromtimestamp(split).isoformat(),'last':datetime.fromtimestamp(pairs[-1][0]).isoformat(),'intercept_C':alpha,'slope':beta,'blocked_test_MAE_C':sum(abs(e) for e in errors)/len(errors),'blocked_test_bias_C':sum(errors)/len(errors),'blocked_test_max_absolute_C':max(abs(e) for e in errors),'blocked_test_RMSE_C':math.sqrt(sum(e*e for e in errors)/len(errors)),'future_transfer_assumption':'Same relationship after weather-tower fall is untested; all-model reconstructed forcing must share policy and uncertainty.'}
Path(__file__).with_name('TOWER3_FEATURE_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
