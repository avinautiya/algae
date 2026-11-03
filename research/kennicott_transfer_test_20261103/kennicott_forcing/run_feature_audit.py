"""Build source-only nominal-clock boundary traces, without any melt prediction."""
import json,pathlib,hashlib
from .boundary import build_boundary
P=pathlib.Path(__file__).parent
PAIRS=[('KDT_2020b','2020-07-25T20:00:00Z','2020-09-06T20:00:00Z'),('KDT_2023a','2023-06-08T20:00:00Z','2023-09-02T20:00:00Z'),('KDT_2023b','2023-06-08T20:00:00Z','2023-09-02T20:00:00Z'),('KDT_2023c','2023-06-08T20:00:00Z','2023-09-02T20:00:00Z')]
def main():
 d=P/'nominal_boundaries';d.mkdir(exist_ok=True);out=[]
 for prof,start,end in PAIRS:
  for air in ['Debris_T2','Debris_T3']:
   b=build_boundary(prof,start,end,air_station=air)
   path=d/(prof+'__'+air+'.json');path.write_text(json.dumps(b,separators=(',',':'))+'\n');row={'file':str(path.relative_to(P)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),**b['metadata']};out.append(row);print(prof,air,'proxyhours',row['proxy_hours'],'fit',row['fit_intercept_C'],row['fit_slope'],'RMSE',row['blocked_RMSE_C'])
 (P/'FEATURE_ONLY_BOUNDARY_AUDIT.json').write_text(json.dumps({'uses_melt_targets':False,'time_policy':'nominalADTnoon assumption only; not final clock-family freeze','boundary_files':out},indent=2)+'\n')
if __name__=='__main__':main()
