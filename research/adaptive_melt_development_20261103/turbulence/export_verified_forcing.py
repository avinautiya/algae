"""Whitelist meteorological inputs only: no ablation or fitted parameters."""
import csv,datetime as dt,hashlib,json
from pathlib import Path
import numpy as np
from scipy.io import loadmat,whosmat
ROOT=Path(__file__).parent
FIELDS=['Date','DeltaGMT','Lat','Lon','Zbas','Pre','Ta','Tdew','Ws','ea','esat','Latm','SAB1','SAB2','SAD1','SAD2','Pr','t_bef','t_aft']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 report=[]
 for tower in [1,2]:
  name=f'run_data_Piramide_Tower{tower}_26092023_01-Nov-2022-24-Apr-2023.mat';p=ROOT/name
  d=loadmat(p,variable_names=FIELDS);dates=[]
  for v in d['Date'].ravel():
   s=(v%1)*86400;assert abs(s-round(s))<1e-4
   dates.append(dt.datetime.fromordinal(int(v))-dt.timedelta(days=366)+dt.timedelta(seconds=round(s)))
  columns={'SOURCE_MODEL_LOCAL_LABEL':[x.isoformat()for x in dates],'MATLAB_DATENUM':d['Date'].ravel().tolist()}
  for field in ['Pre','Ta','Tdew','Ws','ea','esat','Latm','SAB1','SAB2','SAD1','SAD2','Pr']:
   x=d[field].ravel();assert len(x)==len(dates)and np.isfinite(x).all();columns[field]=x.tolist()
  columns['Pre_Pa']=(d['Pre'].ravel()*100).tolist()
  columns['shortwave_partition_component_sum_W_m2']=(d['SAB1']+d['SAB2']+d['SAD1']+d['SAD2']).ravel().tolist()
  out=ROOT/f'Tower{tower}_MODEL_FORCING.csv'
  with out.open('w')as f:
   w=csv.writer(f);w.writerow(columns);w.writerows(zip(*columns.values()))
  report.append({'tower':tower,'source_member':name,'source_sha256':sha(p),'output':out.name,'output_sha256':sha(out),'rows':len(dates),'metadata_only_scalars':{x:float(d[x].item())for x in ['DeltaGMT','Lat','Lon','Zbas','t_bef','t_aft']},'date_decode':'MATLABdatenum converted with integer-second residual<1e-4; no timezone assigned','unit_provenance':{'Pre':'hPa: sourceAerodynamic_Resistence multipliesPreby100 toPa','Ta':'C: samefunctioninputcomment','ea':'Pa: sourceTCdriverDs=esat-ea commentPa','Latm':'incominglongwaveW/m2 per sourceforcing/operator; TCdriverN=Latm','SAB1_SAB2_SAD1_SAD2':'sourceinput solarpartitions, exportedindividually; sum suppliedascomponent-sum not independentlymeasurednetradiation','Ws':'m/s sourcepaper/profile/native convention; independentlysourceprocessedwind','Pr':'sourceprecipitationinput; unitsnotverifiedhere, do notconsumeascalibratedwaterflux'},'camera_or_optimized_parameter_values_loaded':False})
 (ROOT/'FORCING_EXPORT_MANIFEST.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
