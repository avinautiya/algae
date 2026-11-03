"""Execute a source-bound<=11day thermal block; never reads melt outcomes."""
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
from three_sensor import reconstruct,UnsupportedModel

def run(block,**options):
 if not isinstance(block.get('source_sha256'),str)or len(block['source_sha256'])!=64 or any(c not in'0123456789abcdef'for c in block['source_sha256']):raise ValueError('bound sourceSHA256 required')
 dates=[dt.datetime.fromisoformat(x)for x in block['source_local_labels']]
 if any(x.tzinfo is not None for x in dates):raise ValueError('do not invent timezone for naive source')
 times=[(x-dates[0]).total_seconds()for x in dates]
 if not 0<times[-1]<=11*86400:raise UnsupportedModel('sourceblockmustbepositiveand<=11days; no fabricated14daytest')
 r=reconstruct(block['depths_m'],times,block['temperature_C'],**options)
 for trace in r['normalized_heat_trace']:
  trace['start_source_local_label']=(dates[0]+dt.timedelta(seconds=trace['start_seconds'])).isoformat();trace['end_source_local_label']=(dates[0]+dt.timedelta(seconds=trace['end_seconds'])).isoformat()
 h=block['depths_m'][-1]-block['depths_m'][-2];T=block['temperature_C']
 conventional=[{'start_seconds':times[j],'end_seconds':times[j+1],'signed_deep_gradient_K_s_m':(times[j+1]-times[j])*((T[j][-2]-T[j][-1])+(T[j+1][-2]-T[j+1][-1]))/(2*h),'duration_seconds':times[j+1]-times[j]}for j in range(len(times)-1)]
 r.update({'source_sha256':block['source_sha256'],'source_block_id':block.get('source_block_id'),'source_target_values_read':False,'baseline_dose_trace':conventional,'baseline_notes':'Conventional gradient uses same native deepest span; constant rate usesduration. Fixed-conductivityflux is externally supplied k timesgradient. All calibrated comparisons requirethe sameindependentTRAINtarget andexactnativeoperators.'})
 return r

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--ice-depth-m',required=True,type=float);p.add_argument('--basal-temperature-C',required=True,type=float);p.add_argument('--balance-window-hours',type=float,default=6);p.add_argument('--minimum-singular-ratio',type=float,default=.001);p.add_argument('--sensor-absolute-bound-C',type=float);a=p.parse_args()
 raw=Path(a.input).read_bytes();b=json.loads(raw);r=run(b,ice_depth_m=a.ice_depth_m,basal_temperature_C=a.basal_temperature_C,balance_window_seconds=a.balance_window_hours*3600,minimum_singular_ratio=a.minimum_singular_ratio,sensor_absolute_bound_C=a.sensor_absolute_bound_C)
 r['input_artifact_sha256']=hashlib.sha256(raw).hexdigest()
 with Path(a.output).open('x')as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps({'supported':r['fit']['supported_point_estimate'],'normalized_total_interface_heat_K_m':r['normalized_total_interface_heat_K_m'],'negative_heat_steps':r['negative_interface_heat_steps'],'melt_validation_claim':False},indent=2))
