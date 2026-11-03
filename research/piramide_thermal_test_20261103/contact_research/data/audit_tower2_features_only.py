"""Tower2 source eligibility + TRAIN feature audit. No camera values emitted."""
from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd
from matio import load_from_mat
BASE=Path(__file__).parents[2];OUT=Path(__file__).parent
sys.path.insert(0,str(BASE/'data'))
from read_verified_process import read_process
CAM=BASE/'data/camera';source=CAM/'Ablation_Tower2.mat'
manifest=json.loads((CAM/'SEALED_SOURCE_MANIFEST.json').read_text());source_sha=next(x['sha256'] for x in manifest['members'] if x['file']==source.name)
assert hashlib.sha256(source.read_bytes()).hexdigest()==source_sha
obj=load_from_mat(source,raw_data=True)['tracks'].properties['any']
names=[str(x.item()) for x in obj['varNames'][0,0].ravel()];assert names==['distance','ablation']
epoch=np.asarray(obj['rowTimes'][0,0].properties['data']).ravel();clock=pd.DatetimeIndex(epoch.astype('int64').astype('datetime64[ms]'))
finite=np.isfinite(np.asarray(obj['data'][0,0].ravel()[0]).ravel())
nominals=['2022-12-02T03:00','2022-12-16T03:00','2022-12-30T03:00','2023-01-13T03:00','2023-01-27T03:00','2023-02-10T03:00']
bindings=[]
for n in nominals:
 t=pd.Timestamp(n);ids=np.flatnonzero(np.abs((clock-t).total_seconds())<=43200)
 eligible=[i for i in ids if finite[i]]
 chosen=min(eligible,key=lambda i:(abs((clock[i]-t).total_seconds()),clock[i])) if eligible else None
 exact=clock.get_indexer([t])[0]
 bindings.append({'nominal':t.isoformat(),'selected':clock[chosen].isoformat() if chosen is not None else None,
                  'nominal_isfinite':bool(finite[exact]) if exact>=0 else None,
                  'nearby_time_finiteness':[{'label':clock[i].isoformat(),'isfinite':bool(finite[i])} for i in ids]})
thermal=read_process(BASE/'data/decoded','debtemp','Tower2');air=read_process(BASE/'data/decoded','trh','Tower2');moist=read_process(BASE/'data/decoded','sm','Tower2')

def cover(frame,cols,a,b):
 valid=frame[cols].notna().all(axis=1)&np.isfinite(frame[cols]).all(axis=1);idx=frame.index[valid];s,e=pd.Timestamp(a),pd.Timestamp(b)
 left=idx[idx<=s];right=idx[idx>=e];inside=idx[(idx>=s)&(idx<=e)]
 bounded=bool(len(left) and len(right));labels=sorted(set(list(inside)+([left[-1]] if len(left) else [])+([right[0]] if len(right) else [])))
 gaps=[(y-x).total_seconds()/3600 for x,y in zip(labels[:-1],labels[1:]) if y>s and x<e]
 gap=max(gaps,default=0.)
 if not bounded:gap=max(gap,max(0,(idx[0]-s).total_seconds()/3600),max(0,(e-idx[-1]).total_seconds()/3600))
 return {'native_finite_points_inside':len(inside),'source_bounds_cover':bounded,'maximum_gap_hours':gap,'pass_no_extrapolation_max18h':bounded and gap<=18,'first_native_label':idx[0].isoformat(),'last_native_label':idx[-1].isoformat()}

intervals=[]
for j,(a,b) in enumerate(zip(bindings[:-1],bindings[1:])):
 s,e=a['selected'],b['selected'];row={'role':'train' if j==0 else 'heldout','start':s,'end':e}
 if s is None or e is None:row.update({'source_eligible':False,'reason':'No finite camera endpoint withinfixed12h'})
 else:
  row['thermal']=cover(thermal,list(thermal.columns),s,e);row['air_2m']=cover(air,['TA_3'],s,e)
  row['source_eligible']=row['thermal']['pass_no_extrapolation_max18h'] and row['air_2m']['pass_no_extrapolation_max18h']
 intervals.append(row)
s,e=intervals[0]['start'],intervals[0]['end'];assert s and e
sm=moist.loc[(moist.index>=s)&(moist.index<e)];t=thermal.loc[(thermal.index>=s)&(thermal.index<e)];phi=.44
qc={}
for c in sm:
 good=sm[c].between(0,phi)&sm[c].notna();qc[c]={'valid':int(good.sum()),'rows':len(sm),'negative':int((sm[c]<0).sum()),'above_phi':int((sm[c]>phi).sum()),'valid_mean':float(sm.loc[good,c].mean()) if good.any() else None}
complete=((sm>=0)&(sm<=phi)&sm.notna()).all(axis=1);theta=float(sm.loc[complete].mean(axis=1).mean())
rho=2657.3133333333335;cp=816.1266666666667
C=lambda w:(1-phi)*rho*cp+4180000*w+1206*(phi-w)
# Exploratory quasi-steady bottom extrapolation, TRAIN temperatures ONLY.
# Predictor and response share DT5 error; this is not an identifiable physical R estimate.
x=(t.DT4-t.DT5)/.11;y=t.DT5;ell=float(np.dot(x,y)/np.dot(x,x))
out={'scope':'Tower2 feature/time/camera-finiteness ONLY; no camera distance or ablation values read for comparison/scoring',
 'camera_source_sha256':source_sha,'clock':'source_local_naive_no_utc_claim','selection_rule':'For each fixednominal anchor nearestFINITE nativecamera distance timestamp within±12h,tieearlier,nointerpolation; nooutcomefit',
 'endpoint_bindings':bindings,'intervals':intervals,'endpoint_collision':len({b['selected'] for b in bindings if b['selected']})!=sum(b['selected'] is not None for b in bindings),
 'thermal_native_step_seconds':1800,'thermal_native_uniform_step_verified':bool(np.all(np.diff(thermal.index.values).astype('timedelta64[s]').astype(int)==1800)),
 'source_geometry':{'thickness_m':.45,'thermistor_depth_m':[0,.11,.22,.33,.44],'depth_and_thickness_error_m':.02,'bottom_nominal_unsensed_gap_m':.01,'moisture_depth_m':[.19,.29,.39],'tower_fall':'2023-02-12','source':'https://doi.org/10.1029/2025JF008360'},
 'training_features':{'start':s,'end':e,'rows':len(t),'moisture_qc':qc,'all_sm_valid_rows':int(complete.sum()),'equal_probe_theta_proxy_not_volume_mean':theta,
 'conditional_C_J_m3K':{'dry':C(0),'training_proxy':C(theta),'water_saturated':C(phi)},'rho_rock_kg_m3':rho,'cp_rock_J_kgK':cp,'porosity_assumed_not_measured':phi,
 'temperature_min':{c:float(t[c].min()) for c in t},'temperature_max':{c:float(t[c].max()) for c in t},'bottom_above_0p4C_count':int((t.DT5>.4).sum()),
 'quasisteady_effective_extrapolation_length_m':ell,'qs_formula':'OLS T5 = ell*(T4-T5)/.11; temperature-only exploratory diagnosis, sharedsensorerror/storageconfounded, NOTintrinsicRfit',
 'identifiability_warning':'ell combines unsensed basal thickness, kRcontact and transient storage; independent intrinsic contact resistance is not identified by this statistic'},
 'camera_outcome_values_exposed':False,'camera_model_score_computed':False}
(OUT/'TOWER2_FEATURE_ELIGIBILITY.json').write_text(json.dumps(out,indent=2)+'\n')
print('Tower2 feature/finiteness-only audit saved.')
for row in intervals:print(row['role'],row['start'],row['end'],row['source_eligible'])
print('Training capacity scenes',out['training_features']['conditional_C_J_m3K'])
