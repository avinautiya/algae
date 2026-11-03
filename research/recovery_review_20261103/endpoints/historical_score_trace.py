"""Declared-window, piecewise-linear engineering endpoint scorer."""
import argparse,hashlib,json,math
from pathlib import Path

def finite(x):
 try:return type(x) in (int,float) and math.isfinite(x)
 except OverflowError:return False

def validate(trace,spec):
 required={'pre_start_s','pulse_start_s','pulse_end_s','post_end_s','contact_threshold','maximum_sampling_gap_s'}
 if not isinstance(spec,dict) or set(spec)!=required:raise ValueError('exact declared-window specification required')
 if not all(finite(x) for x in spec.values()):raise ValueError('finite scalar specification required')
 if not(spec['pre_start_s']<spec['pulse_start_s']<spec['pulse_end_s']<spec['post_end_s']):raise ValueError('strictly ordered pre/pulse/post windows required')
 if not 0<spec['contact_threshold']<1 or spec['maximum_sampling_gap_s']<=0:raise ValueError('invalid threshold or gap limit')
 if not isinstance(trace,list) or len(trace)<2:raise ValueError('at least two trace records required')
 for row in trace:
  if not isinstance(row,dict) or set(row)!={'time_s','toe_fraction','delivered_heat_W'}:raise ValueError('trace schema mismatch')
  if not all(finite(v) for v in row.values()) or not 0<=row['toe_fraction']<=1:raise ValueError('invalid trace values')
 for a,b in zip(trace,trace[1:]):
  dt=b['time_s']-a['time_s']
  if dt<=0 or dt>spec['maximum_sampling_gap_s']:raise ValueError('nonmonotone time or unsupported sampling gap')
 if trace[0]['time_s']>spec['pre_start_s'] or trace[-1]['time_s']<spec['post_end_s']:raise ValueError('trace does not cover declared windows')

def interpolate(trace,t):
 for a,b in zip(trace,trace[1:]):
  if a['time_s']<=t<=b['time_s']:
   f=(t-a['time_s'])/(b['time_s']-a['time_s'])
   return dict(time_s=t,**{k:a[k]+f*(b[k]-a[k]) for k in ['toe_fraction','delivered_heat_W']})
 raise ValueError('time outside trace')

def window(trace,start,end,threshold):
 rows=[interpolate(trace,start)]+[dict(x) for x in trace if start<x['time_s']<end]+[interpolate(trace,end)]
 # Insert exact linear contact crossings before measuring durations/integrals.
 split=[]
 for a,b in zip(rows,rows[1:]):
  split.append(a);fa=a['toe_fraction']-threshold;fb=b['toe_fraction']-threshold
  if fa*fb<0:
   t=a['time_s']+(b['time_s']-a['time_s'])*(-fa)/(fb-fa);split.append(interpolate(trace,t))
 split.append(rows[-1]);wet=0.;heat=0.;share=0.;events=[]
 for a,b in zip(split,split[1:]):
  dt=b['time_s']-a['time_s'];mid=(a['toe_fraction']+b['toe_fraction'])/2
  if mid>=threshold:wet+=dt
  heat+=dt*(a['delivered_heat_W']+b['delivered_heat_W'])/2
  share+=dt*(a['toe_fraction']+b['toe_fraction'])/2
 for a,b in zip(rows,rows[1:]):
  fa=a['toe_fraction']-threshold;fb=b['toe_fraction']-threshold
  if fa>=0 and fb<0 or fa<0 and fb>=0:
   t=a['time_s']+(b['time_s']-a['time_s'])*(-fa)/(fb-fa)
   events.append(dict(time_s=t,event='loss' if fb<0 else 'reconnection'))
 return dict(wet_seconds=wet,dry_seconds=end-start-wet,signed_delivered_heat_J=heat,time_mean_toe_fraction=share/(end-start),events=events)

def score(trace,spec):
 validate(trace,spec);th=spec['contact_threshold'];start=spec['pre_start_s'];finish=spec['post_end_s'];post=spec['pulse_end_s']
 all_=window(trace,start,finish,th);pw=window(trace,post,finish,th)
 losses=[e for e in all_['events'] if e['event']=='loss'];initially_wet=interpolate(trace,start)['toe_fraction']>=th
 genuine_loss=bool(losses);first_loss=losses[0]['time_s'] if losses else None
 reconnects=[e for e in all_['events'] if e['event']=='reconnection' and first_loss is not None and e['time_s']>first_loss]
 reloss=bool(reconnects and any(e['time_s']>reconnects[0]['time_s'] for e in losses))
 sustained=pw['dry_seconds']<=1e-12*max(1.,finish-post)
 if not genuine_loss and initially_wet:label='maintained_contact'
 elif not genuine_loss:label='initially_disconnected_recovered' if sustained else 'initially_disconnected'
 elif sustained:label='sustained_postpulse_recovery'
 elif reloss:label='reconnected_then_relost'
 elif reconnects:label='temporary_or_incomplete_reconnection'
 else:label='lost_without_reconnection'
 # Hash sampled prefix only: equality establishes duplicate histories, not independence.
 cutoff=first_loss if first_loss is not None else spec['pulse_start_s']
 prefix=[dict(r,time_s=r['time_s']-start) for r in trace if start<=r['time_s']<cutoff]
 r=interpolate(trace,cutoff);r['time_s']-=start;prefix.append(r)
 fingerprint=hashlib.sha256(json.dumps(prefix,sort_keys=True,allow_nan=False,separators=(',',':')).encode()).hexdigest()
 return dict(classification=label,initially_wet=initially_wet,genuine_loss=genuine_loss,first_loss_s=first_loss,reconnected_after_loss=bool(reconnects),reloss=reloss,sustained_over_entire_postpulse_window=sustained,postpulse=pw,whole_window=all_,sampled_preloss_history_sha256=fingerprint,specification_sha256=hashlib.sha256(json.dumps(spec,sort_keys=True,separators=(',',':')).encode()).hexdigest(),interpretation='Declared engineering endpoints under linear interpolation; not statistical confidence, independent histories, field validation or discovery certification.')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('trace',type=Path);p.add_argument('spec',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 if a.out.exists():p.error('output exists; choose a versioned new file')
 def reject(x):raise ValueError('nonfinite JSON constant '+x)
 def unique(pairs):
  d={}
  for k,v in pairs:
   if k in d:raise ValueError('duplicate JSON key '+k)
   d[k]=v
  return d
 r=score(json.loads(a.trace.read_text(),parse_constant=reject,object_pairs_hook=unique),json.loads(a.spec.read_text(),parse_constant=reject,object_pairs_hook=unique));a.out.parent.mkdir(parents=True,exist_ok=True)
 with a.out.open('x') as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
if __name__=='__main__':main()
