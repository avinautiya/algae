"""Typed process features for preregistration; no camera inputs or target scores."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
P=Path(__file__).parent;D=Path('/private/tmp/glacier-round2-20261103/data/decoded');sys.path.insert(0,str(D.parent))
from read_verified_process import read_process
A0=pd.Timestamp('2022-12-02T03:00');A1=pd.Timestamp('2023-02-10T03:00')
manifest={'scope':'Canonical observedprocessfeatures only; noTower2cameraoutcomes','clock':'source_local_naive_no_utc_claim',
 'air_operator':'TA_3labelvalues linearlyinterpolated on nativeSM3halfhourlabels; requirefinitebrackets,max18h; noextrapolation; source-row averaging/labelphase unresolved',
 'soil_T_operator':{'Tower1':'DT5at.50m exactlynominalSM3depth','Tower2':'linearDT4(.33m)toDT5(.44m)atSM3(.39m); derivedtemperatureQC, notco-locatedmeasurement'},
 'temperature_accuracy_C':.4,'manufacturer_moisture_typical_absolute_accuracy_m3_m3':.031,'localcalibrationcertificate':None,'sites':{}}
for j in [1,2]:
 theta=read_process(D,'sm',f'Tower{j}')['SM3'];temp=read_process(D,'debtemp',f'Tower{j}');air=read_process(D,'trh',f'Tower{j}')['TA_3']
 theta=theta.loc[(theta.index>=A0)&(theta.index<=A1)];temp=temp.loc[theta.index]
 source=air.loc[(air.index>=A0-pd.Timedelta(hours=18))&(air.index<=A1+pd.Timedelta(hours=18))]
 finite=source.notna()&np.isfinite(source);source=source.loc[finite]
 sn=source.index.as_unit('ns').asi8;tn=theta.index.as_unit('ns').asi8
 right=np.searchsorted(sn,tn,side='left');left=right-1
 exact=(right<len(sn))&(sn[np.minimum(right,len(sn)-1)]==tn);left[exact]=right[exact]
 if (left<0).any() or (right>=len(sn)).any():raise ValueError('Sourceboundsnotcovercanonicalgrid')
 gap=(sn[right]-sn[left])/3.6e12
 if (gap>18).any():raise ValueError('Airinterpolationgap>18h')
 vals=np.interp(tn.astype(float),sn.astype(float),source.to_numpy())
 soil=temp.DT5 if j==1 else temp.DT4+(temp.DT5-temp.DT4)*(.39-.33)/(.44-.33)
 frame=pd.DataFrame({'TA_3_C_interpolated':vals,'SM3_raw_m3_m3':theta.to_numpy(),'soilT_at_SM3_C_reference':soil.to_numpy(),
  'soilT_reference_below0':(soil<0).to_numpy(),'soilT_reference_above50':(soil>50).to_numpy(),'soilT_reference_definitely_within_0to50_given_pminus0p4':((soil>=.4)&(soil<=49.6)).to_numpy(),
  'SM3_negative':(theta<0).to_numpy(),'SM3_above_assumed_phi0p44':(theta>.44).to_numpy(),'air_bracket_gap_hours':gap,'air_exact_native_label':exact,
  'air_left_source_label':source.index[left].astype(str),'air_right_source_label':source.index[right].astype(str)},index=theta.index)
 if not np.isfinite(frame[['TA_3_C_interpolated','SM3_raw_m3_m3','soilT_at_SM3_C_reference']]).all().all():raise ValueError('Nonfinitecanonicalfeatures')
 file=P/f'Tower{j}_CANONICAL_FEATURES.csv';frame.to_csv(file,index_label='SOURCE_LOCAL_LABEL',date_format='%Y-%m-%dT%H:%M:%S')
 af=P/f'Tower{j}_NATIVE_AIR_TA3.csv';source.to_csv(af,index_label='SOURCE_LOCAL_LABEL',date_format='%Y-%m-%dT%H:%M:%S')
 train=frame.loc[(frame.index>=A0)&(frame.index<pd.Timestamp('2022-12-16T03:00'))]
 anchors=['2022-12-02T03:00:00','2022-12-16T03:00:00', '2022-12-30T00:00:00' if j==1 else '2022-12-30T03:00:00','2023-01-13T03:00:00','2023-01-27T03:00:00','2023-02-10T03:00:00']
 site={'rows':len(frame),'first_label':frame.index[0].isoformat(),'last_label':frame.index[-1].isoformat(),'training_SM3_mean_halfopen':float(train.SM3_raw_m3_m3.mean()),'training_rows':len(train),
   'frozen_site_specific_anchors':anchors,'target_role':'AlreadyinspectedT1development' if j==1 else 'CameraT2unopenedcandidateevaluation',
   'soilT_QC_counts':{'nominal_below0':int(frame.soilT_reference_below0.sum()),'nominal_above50':int(frame.soilT_reference_above50.sum()),'not_guaranteed_0to50_given_pminus0p4':int((~frame.soilT_reference_definitely_within_0to50_given_pminus0p4).sum())},
   'source_sha256':{f'{f}_Tower{j}.csv':hashlib.sha256((D/f'{f}_Tower{j}.csv').read_bytes()).hexdigest() for f in ['sm','debtemp','trh']},
   'canonical_feature_file':file.name,'canonical_feature_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'native_air_file':af.name,'native_air_sha256':hashlib.sha256(af.read_bytes()).hexdigest(),
   'process_source_manifest':'/private/tmp/glacier-round2-20261103/data/decoded/DECODE_VERIFICATION.json'}
 manifest['sites'][f'Tower{j}']=site
manifest['builder_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(P/'CANONICAL_FEATURE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Two canonical featurestreams built; no camera outcomes accessed.')
for k,v in manifest['sites'].items():print(k,v['rows'],v['training_SM3_mean_halfopen'],v['soilT_QC_counts'])
