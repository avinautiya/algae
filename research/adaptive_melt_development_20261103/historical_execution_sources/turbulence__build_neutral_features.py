"""Build weather/profile-only features over already-released development windows."""
import bisect,csv,datetime as dt,hashlib,json,math
from pathlib import Path
from neutral_profile import neutral_flux
ROOT=Path(__file__).parent;OLD=Path('/private/tmp/glacier-round2-20261103/data/decoded')
def readcsv(p):
 with Path(p).open()as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 reports=[]
 for tower in [1,2]:
  trh=readcsv(OLD/f'trh_Tower{tower}.csv');winds={r['SOURCE_LOCAL_LABEL']:r for r in readcsv(OLD/f'wind_Tower{tower}.csv')};pr=readcsv(ROOT/f'Tower{tower}_MODEL_FORCING.csv');pt=[dt.datetime.fromisoformat(r['SOURCE_MODEL_LOCAL_LABEL'])for r in pr];pv=[float(r['Pre_Pa'])for r in pr]
  output=[]
  for row in trh:
   t=dt.datetime.fromisoformat(row['SOURCE_LOCAL_LABEL'])
   if not dt.datetime(2022,12,2,3)<=t<=dt.datetime(2023,2,10,3):continue
   wind=winds.get(row['SOURCE_LOCAL_LABEL'])
   if wind is None:raise ValueError('no native wind match')
   idx=bisect.bisect_right(pt,t)-1
   if idx<0 or idx+1>=len(pt):raise ValueError('no pressure extrapolation')
   a=(t-pt[idx]).total_seconds()/(pt[idx+1]-pt[idx]).total_seconds();pressure=pv[idx]+a*(pv[idx+1]-pv[idx])
   try:r=neutral_flux([float(wind[f'L{i}_FF'])for i in range(1,5)],[float(row[f'TA_{i}'])for i in range(1,5)],[float(row[f'RH_{i}'])for i in range(1,5)],pressure)
   except ValueError as e:r={'status':'unsupported_physical_source_input','source_reason':str(e)}
   output.append({'SOURCE_LOCAL_LABEL':t.isoformat(),**r})
  names=['SOURCE_LOCAL_LABEL']+sorted(set().union(*(r.keys()for r in output)) - {'SOURCE_LOCAL_LABEL'})
  out=ROOT/f'Tower{tower}_NEUTRAL_PROFILE_FEATURES.csv'
  with out.open('w')as f:w=csv.DictWriter(f,fieldnames=names);w.writeheader();w.writerows(output)
  valid=[r for r in output if r.get('H_up_W_m2')is not None]
  reports.append({'tower':tower,'feature_csv':out.name,'sha256':sha(out),'rows':len(output),'supported_neutral_proxy_rows':len(valid),'refused_rows':len(output)-len(valid),'positive_wind_slope_robust_given_accuracy':sum(r.get('wind_slope_sign_resolved_given_stated_accuracy',False)for r in output),'temperature_slope_sign_resolved_given_accuracy':sum(r.get('T_slope_sign_resolved_given_stated_accuracy',False)for r in output),'rows_abs_gradientRi_le_point03':sum(abs(r['gradient_Ri_at1m_approx'])<=.03 for r in valid),'H_range_W_m2':[min(r['H_up_W_m2']for r in valid),max(r['H_up_W_m2']for r in valid)],'LE_range_W_m2':[min(r['LE_up_W_m2']for r in valid),max(r['LE_up_W_m2']for r in valid)],'source_hashes':{str(p):sha(p)for p in [OLD/f'trh_Tower{tower}.csv',OLD/f'wind_Tower{tower}.csv',ROOT/f'Tower{tower}_MODEL_FORCING.csv']},'clock_operator':'sourceauthornaivelabelbinding; pressurelinearlyinterpolatedbetweenhourlymodelinputlabels; no target-dependent lag','no_melt_target_access':True,'actual_turbulent_flux_validated':False,'closure':'allfourheightsOLSagainstlnz; neutral,d=0,constantflux; upwardpositive; nonpositivewindslopeexplicitrefusal'})
 (ROOT/'NEUTRAL_FEATURE_MANIFEST.json').write_text(json.dumps(reports,indent=2)+'\n');print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
