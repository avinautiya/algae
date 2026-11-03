"""Process-feature audit, both towers. Tower2 camera outcomes never accessed."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
P=Path(__file__).parent;sys.path.insert(0,'/private/tmp/glacier-round2-20261103/data')
from read_verified_process import read_process
D=Path('/private/tmp/glacier-round2-20261103/data/decoded')
anchors=pd.to_datetime(['2022-12-02T03:00','2022-12-16T03:00','2022-12-30T03:00','2023-01-13T03:00','2023-01-27T03:00','2023-02-10T03:00'])
w=pd.read_csv(P/'PIR22-23_60-kd_ERA5.csv');w.index=pd.to_datetime(w.DATE);w=w.drop(columns='DATE')
out={'scope':'T1development/T2processfeatures only, NO Tower2camera values or scores','source_weather_sha256':hashlib.sha256((P/'PIR22-23_60-kd_ERA5.csv').read_bytes()).hexdigest(),
 'weather_provenance':'Published calibration input, ERA5 in filename; mixed measured/filled/lapse-corrected provenance, not rawcolocatedprecip/snowmeasurements',
 'period_definition':'Common half-open nominalT2calendar windows, allT03; Tower1earlierfieldendDec30T00 differs3h, no comparisonmetricrecomputedhere',
 'manufacturer_accuracy':'S-SMC-M005 typicalabsolute±.031m3/m3 formineralsoil0..50C; not3percentofreading and not guaranteedrock-debrisfieldbound',
 'towers':{}}
for j,d,depths in [(1,.5,[.3,.4,.5]),(2,.45,[.19,.29,.39])]:
 sm=read_process(D,'sm',f'Tower{j}');T=read_process(D,'debtemp',f'Tower{j}')
 sm0=sm.loc[(sm.index>=anchors[0])&(sm.index<anchors[1])].mean()
 rows=[]
 for n,(a,b) in enumerate(zip(anchors[:-1],anchors[1:])):
  s=sm.loc[(sm.index>=a)&(sm.index<b)];t=T.loc[(T.index>=a)&(T.index<b)];r={'role':'train' if n==0 else 'future_feature_only','start':a.isoformat(),'end':b.isoformat(),'rows_sm':len(s),'rows_T':len(t),'probe_stats':{},'temperature_stats':{}}
  for c in s:
   valid=s[c].between(0,.44)&s[c].notna()
   r['probe_stats'][c]={'raw_mean':float(s[c].mean()),'raw_min':float(s[c].min()),'raw_max':float(s[c].max()),'negative':int((s[c]<0).sum()),'above_assumed_phi_0p44':int((s[c]>.44).sum()),'missing':int(s[c].isna().sum()),'physically_valid_mean':float(s.loc[valid,c].mean()) if valid.any() else None,'ratio_to_own_training_mean':float(s[c].mean()/sm0[c]),'change_from_own_training_mean':float(s[c].mean()-sm0[c]),'below_manufacturer_absolute_accuracy':int((s[c]<.031).sum())}
  for c in ['DT3','DT4','DT5']:
   r['temperature_stats'][c]={'min_C':float(t[c].min()),'max_C':float(t[c].max()),'below_0C':int((t[c]<0).sum()),'below_minus0p4C':int((t[c]<-.4).sum()),'above_50C':int((t[c]>50).sum())}
  ww=w.loc[(w.index>=a)&(w.index<b)]
  r['shared_weather_file']={'rows':len(ww),'PP_sum_mm_as_supplied':float(ww.PP.sum()) if ww.PP.notna().any() else None,'PP_missing':int(ww.PP.isna().sum()),'PP_positive_hours':int((ww.PP>0).sum()),'SNOW_missing':int(ww.SNOW.isna().sum()),'SNOW_positive':int((ww.SNOW>0).sum()),'albedo_directly_observed':False,'SWIN_available':int(ww.SWIN.notna().sum()),'SWOUT_available':int(ww.SWOUT.notna().sum()),'source_clock_assumption':'sameunlabeledlocalclock; rowlabelmeanalignmentnotindependentlyresolved'}
  rows.append(r)
 out['towers'][f'Tower{j}']={'debris_thickness_m':d,'probe_depths_m':depths,'nominal_probe_gaps_to_ice_m':[round(d-z,3) for z in depths],'probe_gap_separate_error_box_m':[ [max(0.,round(d-z-.04,3)),round(d-z+.04,3)] for z in depths], 'depth_thickness_uncertainty_each_m':.02,'source_file_sha256':{x:hashlib.sha256((D/f'{x}_Tower{j}.csv').read_bytes()).hexdigest() for x in ['sm','debtemp']},'periods':rows}
(P/'MOISTURE_TRANSFER_FEATURE_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
for tower,v in out['towers'].items():print(tower,'SM3means',[r['probe_stats']['SM3']['raw_mean'] for r in v['periods']])
print('No Tower2 camera data accessed.')
