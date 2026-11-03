"""Reproduce a declared TRAIN sensible-storage proxy; no camera files read."""
import argparse,csv,json
from pathlib import Path

def calculate(process):
 with Path(process).open() as stream:rows={r['SOURCE_LOCAL_LABEL']:r for r in csv.DictReader(stream)}
 depths=[.02,.125,.25,.375,.5,.52];C=1855487.1033697424;duration=14*86400;endpoints={}
 for stamp in ['2022-12-02T03:00:00','2022-12-16T03:00:00']:
  values=[float(rows[stamp][f'DT{i}']) for i in range(1,6)]+[0.]
  integral=sum((b-a)*(ta+tb)/2 for a,b,ta,tb in zip(depths,depths[1:],values,values[1:]))
  endpoints[stamp]=dict(assumed_profile_integral_K_m=integral)
 a,b=endpoints.values();change=C*(b['assumed_profile_integral_K_m']-a['assumed_profile_integral_K_m'])
 return dict(scope='TRAIN endpoint sensible-storage proxy only; piecewise-linear sensor profile and0C base assumptions, not a measured total energy bound',endpoints=endpoints,bulk_capacity_reference_J_m3K=C,storage_change_J_m2=change,mean_storage_change_W_m2=change/duration,conditional_plusminus04C_sensor_induced_mean_storage_radius_W_m2=2*.4*C*.5/duration,ice_storage_and_latent_pore_terms_measured=False,Tower2_camera_read=False)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('process');args=p.parse_args();print(json.dumps(calculate(args.process),indent=2))
