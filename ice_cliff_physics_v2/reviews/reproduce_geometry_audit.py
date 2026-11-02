"""Reproduce the declared synthetic height audit; no field targets or tuning."""
import argparse,copy,hashlib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--package',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 sys.path.insert(0,str(a.package));from model import forecast
 base=json.loads((ROOT/'geometry_counterexample_config.json').read_text());rows=[]
 for count in [4,16,64,128,256]:
  c=copy.deepcopy(base);split=count//4;expanded=[]
  for b in c['upgrades']['initial_face_bands']:
   z=[b['z_low_m']+(b['z_high_m']-b['z_low_m'])*i/split for i in range(split)]+[b['z_high_m']]
   expanded.extend(dict(z_low_m=lo,z_high_m=hi,front_x_m=b['front_x_m']) for lo,hi in zip(z,z[1:]))
  c['upgrades']['face_profile_bands']=count;c['upgrades']['initial_face_bands']=expanded;c['step_s']=300.;c['periods'][0]['duration_s']=21600.
  result=forecast(c);last=result['series'][-1]
  rows.append(dict(n_bands=count,totals=result['totals'],checks=result['checks'],basal_front_m=last['geometry']['face_x_m'],floor_z_m=last['geometry']['floor_z_m'],max_undercut_m=last['undercut_diagnostic']['maximum_horizontal_undercut_m']))
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(dict(scope='synthetic integrated height refinement; no field accuracy',code_sha256=hashlib.sha256((a.package/'model.py').read_bytes()).hexdigest(),input_sha256=hashlib.sha256((ROOT/'geometry_counterexample_config.json').read_bytes()).hexdigest(),rows=rows),indent=2)+'\n')
 print('Wrote',a.out)
if __name__=='__main__':main()
