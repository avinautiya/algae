"""Prospective endpoint taxonomy v2; historical v1 scores stay unchanged.

Uses the same explicitly piecewise-linear observation operator as v1. This
does not certify event timing at discontinuous forcing jumps; those require
one-sided observations or an independently defined jump-aware operator.
"""
import hashlib,json,re
if __package__:
 from . import historical_score_trace as v1
else:
 import historical_score_trace as v1

def score(trace,spec,*,trace_interpretation):
 if trace_interpretation!='piecewise_linear':raise ValueError('explicit piecewise-linear observation contract required')
 v1.validate(trace,spec)
 start=spec['pre_start_s'];pstart=spec['pulse_start_s'];pend=spec['pulse_end_s'];end=spec['post_end_s'];th=spec['contact_threshold']
 whole=v1.window(trace,start,end,th);pulse=v1.window(trace,pstart,pend,th);post=v1.window(trace,pend,end,th)
 initial=v1.interpolate(trace,start)['toe_fraction']>=th;endwet=v1.interpolate(trace,end)['toe_fraction']>=th
 losses=[e['time_s'] for e in whole['events'] if e['event']=='loss'];rises=[e['time_s'] for e in whole['events'] if e['event']=='reconnection']
 firstwet=start if initial else (rises[0] if rises else None)
 firstloss=losses[0] if losses else None
 recon=[t for t in rises if firstloss is not None and t>firstloss]
 reloss=bool(recon and any(t>recon[0] for t in losses))
 sustained=post['dry_seconds']<=1e-12*max(1,end-pend)
 everwet=whole['wet_seconds']>0
 if not initial:
  if not everwet:label='initially_disconnected_never_wet'
  elif sustained:label='initially_disconnected_sustained_postpulse_access'
  elif not endwet:label='initially_disconnected_temporarily_rewetted_then_dry'
  else:label='initially_disconnected_partial_postpulse_access'
 elif not losses:label='maintained_contact'
 elif sustained:label='sustained_postpulse_recovery'
 elif not recon:label='lost_without_reconnection'
 elif reloss:label='reconnected_then_relost'
 else:label='temporary_or_incomplete_reconnection'
 return dict(schema_version=2,classification=label,initially_wet=initial,first_threshold_wetting_s=firstwet,positive_duration_contact_observed=everwet,first_loss_s=firstloss,reconnected_after_loss=bool(recon),reloss_after_reconnection=reloss,wet_at_end=endwet,rewet_during_pulse=pulse['wet_seconds']>0 and not initial,sustained_entire_postpulse_window=sustained,whole_window=whole,pulse=pulse,postpulse=post,
  historical_v1_score=v1.score(trace,spec),specification_sha256=hashlib.sha256(json.dumps(spec,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
  interpretation='Prospective engineering taxonomy. V1 predictions and scores remain authoritative for their frozen evaluation. V2 never converts a historical label mismatch into a success. No field validation or statistical independence is certified.')

def compare_frozen_prediction(*,predicted_classes,historical_class):
 if not isinstance(predicted_classes,list) or not predicted_classes or not all(isinstance(x,str) for x in predicted_classes) or not isinstance(historical_class,str):raise ValueError('exact historical labels required')
 return dict(historical_label_match=historical_class in predicted_classes,retrospective_relabeling_permitted=False)

def sampled_history_digest(historical_score):
 """Use the actual v1 key; missing identities never become null duplicates."""
 value=historical_score.get('sampled_preloss_history_sha256') if isinstance(historical_score,dict) else None
 if not isinstance(value,str) or re.fullmatch('[0-9a-f]{64}',value) is None:raise ValueError('sampled preloss history digest missing; cannot infer independence')
 return value
