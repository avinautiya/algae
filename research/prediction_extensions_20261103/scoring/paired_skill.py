"""Matched-support deterministic daily melt scores; no automatic accuracy claim."""
import json,math

def finite(x):
 try:return type(x) in (float,int) and math.isfinite(x)
 except OverflowError:return False

def score(rows,*,candidate,models,unit,expected_keys):
 """expected_keys is the frozen site/date support, not selected after scoring.

 No missing prediction is dropped. The comparison must be regenerated on a new
 explicitly versioned support if any method cannot predict an eligible record.
 Units/time/observation definitions must be harmonized upstream. Day wins are
 descriptive counts, not independent replications or significance tests.
 """
 if not isinstance(rows,list) or not rows:raise ValueError('nonempty record list required')
 if unit!='mm_water_equivalent_per_day':raise ValueError('explicit common daily melt unit required')
 if not isinstance(models,list) or len(models)<2 or not all(isinstance(m,str) and m.strip() for m in models) or len(set(models))!=len(models) or candidate not in models:raise ValueError('candidate and at least one distinct frozen baseline required')
 if not isinstance(expected_keys,list) or not expected_keys or not all(isinstance(k,str) and k.strip() for k in expected_keys) or len(set(expected_keys))!=len(expected_keys):raise ValueError('unique nonempty string support keys required')
 seen=set();err={m:[] for m in models};sites={};wins={m:dict(candidate_lower_error=0,tie=0,candidate_higher_error=0) for m in models if m!=candidate}
 for row in rows:
  if not isinstance(row,dict):raise ValueError('record object required')
  key=row.get('key');site=row.get('site')
  if key in seen or key not in expected_keys or not isinstance(site,str) or not site:raise ValueError('duplicate/unexpected/unlabelled support')
  seen.add(key)
  obs=row.get('observed');pred=row.get('predictions')
  if not finite(obs) or not isinstance(pred,dict) or set(pred)!=set(models) or not all(finite(v) for v in pred.values()):raise ValueError('finite complete predictions required; no silent deletion')
  if row.get('unit')!=unit:raise ValueError('mixed units')
  # Signed observations/predictions remain visible. Whether freezing/accumulation
  # belongs to the target is an upstream observation-contract decision.
  for m in models:
   delta=pred[m]-obs
   if not finite(delta):raise ValueError('error arithmetic overflows')
   err[m].append(delta)
  sites[site]=sites.get(site,0)+1
  ce=abs(pred[candidate]-obs)
  for m in wins:
   delta=ce-abs(pred[m]-obs)
   wins[m]['tie' if delta==0 else ('candidate_lower_error' if delta<0 else 'candidate_higher_error')]+=1
 if seen!=set(expected_keys):raise ValueError('incomplete frozen support')
 metrics={m:dict(n=len(e),mae=sum(abs(v)/len(e) for v in e),bias=sum(v/len(e) for v in e),rmse=math.hypot(*e)/math.sqrt(len(e))) for m,e in err.items()}
 if not all(finite(v) for d in metrics.values() for v in d.values()):raise ValueError('score arithmetic overflows')
 equal_site={};site_mae={site:{} for site in sites}
 for m in models:
  per_site=[]
  for site in sites:
   values=[abs(r['predictions'][m]-r['observed']) for r in rows if r['site']==site];value=sum(v/len(values) for v in values);per_site.append(value);site_mae[site][m]=value
  equal_site[m]=sum(v/len(per_site) for v in per_site)
 return dict(unit=unit,metrics=metrics,equal_site_mae=equal_site,site_mae=site_mae,day_win_counts=wins,site_counts=sites,candidate_lower_pooled_mae_than_every_baseline=all(metrics[candidate]['mae']<metrics[m]['mae'] for m in wins),candidate_lower_mae_at_every_site_than_every_baseline=all(site_mae[s][candidate]<site_mae[s][m] for s in sites for m in wins),interpretation='Descriptive matched-record comparison. No independence, blind evaluation, confidence interval, novelty or environmental validity is certified.')
