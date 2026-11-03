"""Thermal CSV only. Separate106-byte melt file is not fetched or read."""
from pathlib import Path
import csv,json,hashlib
import pandas as pd,numpy as np
P=Path(__file__).parent;f=P/'IMJA_TEMPERATURE.csv';r=json.loads((P/'IMJA_RECORD.json').read_text());m=next(x for x in r['files'] if x['key']=='15.03743_rounce2015-ts.csv')
assert m['checksum']=='md5:'+hashlib.md5(f.read_bytes()).hexdigest()
rows=list(csv.reader(f.read_text(encoding='utf-8-sig').splitlines()));groups=rows[1];depths=rows[2]
site=None;names=[];mapping=[]
for j in range(1,len(depths)):
 if groups[j].strip():site=groups[j].strip().replace(' ','')
 depth=float(depths[j])/100;name=f'{site}_z{depth:g}m_C';names.append(name);mapping.append({'name':name,'source_site':site,'source_depth_m':depth,'source_unit':'K','converted_unit':'C'})
body=[row for row in rows[3:] if row[0].strip()]
undated_rows_excluded=len(rows[3:])-len(body)
d=pd.DataFrame(body);assert d.shape[1]==24
idx=pd.DatetimeIndex(pd.to_datetime(d[0],format='%m/%d/%y %H:%M'))
assert idx.is_unique and idx.is_monotonic_increasing
v=d.iloc[:,1:].apply(pd.to_numeric,errors='raise').to_numpy()-273.15
frame=pd.DataFrame(v,index=idx,columns=names);frame.to_csv(P/'IMJA_NATIVE_THERMAL.csv',index_label='SOURCE_NEPAL_LABEL')
report={'source_doi':'10.5281/zenodo.3906947','source_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'md5_verified':True,
 'undated_rows_excluded_without_guessing_dates':undated_rows_excluded,'rows':len(idx),'first_label':idx[0].isoformat(),'last_label':idx[-1].isoformat(),'clock':'CSV saysNEPAL; preserve labels,noUTCconversionneededforrelativeanalysis','mapping':mapping,
 'source_accuracy_C':.3,'step_seconds_counts':{str(k):int(n) for k,n in zip(*np.unique(np.diff(idx.values).astype('timedelta64[s]').astype(int),return_counts=True))},
 'missing_per_column':frame.isna().sum().to_dict(),'IM14_bottom_nominal_z0p24m_min_C':float(frame['Site14_z0.24m_C'].min()),
 'IM14_bottom_nominal_z0p24m_below_minus1_count':int((frame['Site14_z0.24m_C']<-1).sum()),
 'source_Rowan_table1_IM14_thickness_m':.26,'geometry_warning':'Deepestsourceprobe.24mwhileRowanTable1layer.26m; sourcecallsinterfacebutdonotsilentlyerase2cmgap',
 'ice_probe_below_interface_present':False,'independent_bulk_C_present':False,'separate_melt_file_opened':False,
 'validation_scope':'Nativeprofiletemperatureavailable; longintervalmeltmetadataexistsbutoutcomesremainsealed, nosufficientwithin-sitechronologicalmeltcalibrationidentified'}
(P/'IMJA_THERMAL_ELIGIBILITY.json').write_text(json.dumps(report,indent=2)+'\n')
print('Imja thermal verified',len(idx),'rows',idx[0],idx[-1]);print('IM14deepmin',report['IM14_bottom_nominal_z0p24m_min_C'],'cold<-1samples',report['IM14_bottom_nominal_z0p24m_below_minus1_count'])
