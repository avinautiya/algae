"""DEVELOPMENT ONLY: all trials retained, no uninspected targets."""
import csv,json,datetime as dt,math
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from package_paths import round2_root,round3_root
from storage_flux import slab_sensible_J_m2,predict
all_sites=[]
old=round2_root();new=round3_root();L=333700.;C=1.855487e6
T2p=json.load(open(new/'moisture/frozen_tower2_final/PREDICTIONS.json'));T2r=json.load(open(new/'moisture/frozen_tower2_final/RESULTS.json'))
t1p=json.load(open(old/'model/frozen_piramide_v3/melt_comparison/PREDICTIONS.json'));t1r=json.load(open(old/'model/frozen_piramide_v3/melt_comparison/RESULTS.json'));t1b=json.load(open(old/'model/frozen_piramide_future_batch_final/melt_comparison/RESULTS.json'))
for tower in [1,2]:
 name='Tower'+str(tower);manifest=json.load(open(new/'data/moisture/CANONICAL_FEATURE_MANIFEST.json'));epochs=manifest['sites'][name]['frozen_site_specific_anchors'];parsed=[dt.datetime.fromisoformat(x)for x in epochs]
 with open(old/f'data/decoded/debtemp_{name}.csv')as f:rows=list(csv.DictReader(f))
 temp={dt.datetime.fromisoformat(x['SOURCE_LOCAL_LABEL']):x for x in rows};spacing=.125 if tower==1 else .11
 if tower==1:
  y0=t1p['calibration']['observed_ice_ablation_mmwe'];n=t1r['primary_candidate'];obs=[t1r['model_metrics'][n]['observed_total_mmwe']]+[x['observed_mmwe']for x in t1b['model_metrics'][n]['period_errors']]
 else:y0=T2p['local_train_target']['observed_ice_ablation_mmwe'];obs=[x['observed_mmwe']for x in T2r['model_metrics']['anomaly_exp_air']['period_errors']]
 Js=[];rawends=[]
 for a,b in zip(parsed[:-1],parsed[1:]):
  times=sorted(t for t in temp if a<=t<=b)
  assert times[0]==a and times[-1]==b and all((v-u).total_seconds()==1800 for u,v in zip(times[:-1],times[1:]))
  gradients=[(float(temp[t]['DT4'])-float(temp[t]['DT5']))/spacing for t in times]
  Js.append(math.fsum((v-u).total_seconds()*(x+y)/2 for u,v,x,y in zip(times[:-1],times[1:],gradients[:-1],gradients[1:])))
  rawends.append((temp[a],temp[b]))
 output=[]
 for capacity in [1.5e6,C,2.5e6]:
  for gap in [0.,.01,.02,.04]:
   def H(row):return slab_sensible_J_m2(float(row['DT4']),float(row['DT5']),spacing,gap,capacity)
   ds=[H(b)-H(a)for a,b in rawends];model=predict(y0,Js,ds);k=model['calibrated_conductivity_W_mK'];pred=model['interval_net_melt_equivalent_mmwe'][1:];errors=[x-y for x,y in zip(pred,obs)]
   output.append({'C':capacity,'gap':gap,'k':k,'MAE':sum(abs(x)for x in errors)/4,'pred':pred,'storage_changes_mmwe':[d/L for d in ds]})
 baseline=[y0*j/Js[0]for j in Js[1:]]
 all_sites.append({'tower':tower,'TRAIN':y0,'obs':obs,'signed_trapezoid_gradient_baseline_MAE':sum(abs(p-y)for p,y in zip(baseline,obs))/4,'all_trials':output})

result={'scope':'T1/T2 released outcomes DEVELOPMENT ONLY; all24 trials retained','no_Tower3_or_later_T1_outcomes_read':True,'sites':all_sites}
path=Path(__file__).with_name('DEVELOPMENT_STORAGE_RESULTS.json')
with path.open('x')as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print(path)
