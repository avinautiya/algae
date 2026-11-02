"""Fail-closed metadata checker, not a scientific-validation certificate.

Version2 intentionally requires typed value/evidence records instead of the
old truthiness test. A complete record still cannot authorize an environmental
accuracy claim. No source is fetched, authenticated or independently reviewed
by this checker. Exit0 means metadata complete, never validated science.
"""
import argparse
import datetime
import hashlib
import ipaddress
import json
import math
import re
from pathlib import Path
from urllib.parse import urlparse

RULES={
 'geometry':dict(floor_elevation_survey_dates='dates',floor_survey_sigma_m='nonnegative',cross_section_width_source='text',width_measured_at_sensor_section='true',sensor_clearance_m_or_on_bed_verified='clearance',depth_basis='text',section_to_station_binding='text'),
 'temperature':dict(bulk_water_temperature_or_validated_floor_to_bulk_relation='true',sensor_calibration_offset_applied_once='true',clock_convention_resolved='true',observed_support_fraction='fraction',reconstructed_support_labelled='true'),
 'discharge':dict(discharge_basis='text',discharge_at_target_section_or_conservation_justified='true',discharge_uncertainty='discharge_sigma'),
 'contact':dict(wetted_contact_observed_or_modelled='contact',missing_record_not_treated_as_zero='true',ice_temperature_known_or_bounded='ice_bounds'),
 'measurement':dict(target_definition_matches_observation_operator='true',target_uncertainty_is_source_sigma='true',independent_of_calibration_data='true',heat_transfer_closure_not_fitted_to_target_event='true'),
 'numerics':dict(solver_version='text',source_sha256='sha256',restart_fingerprint_verified='true',grid_contact_timing='grid',hydraulic_convergence_per_interval='flags',mass_relative_residual='nonnegative',energy_relative_residual='nonnegative',acceptance_limits_predeclared='limits')}

def number(x,minimum=0,maximum=math.inf):
 if type(x) not in (int,float):return False
 try:return math.isfinite(x) and minimum<=x<=maximum
 except (OverflowError,ValueError,TypeError):return False

def text(x):
 return isinstance(x,str) and bool(x.strip()) and x.strip().casefold() not in {'unknown','false','true','none','null','n/a','na','todo','unverified'}

def sha(x):
 return isinstance(x,str) and re.fullmatch('[0-9a-f]{64}',x) is not None and len(set(x))>1

def exact(x,keys):
 return isinstance(x,dict) and set(x)==set(keys)

def valid_value(kind,v):
 if kind=='true':return type(v) is bool and v
 if kind=='text':return text(v)
 if kind=='sha256':return sha(v)
 if kind=='nonnegative':return number(v)
 if kind=='fraction':return number(v,maximum=1)
 if kind=='contact':return v in ('observed','modelled')
 if kind=='flags':return isinstance(v,list) and bool(v) and all(type(x) is bool and x for x in v)
 if kind=='dates':
  if not isinstance(v,list) or not v:return False
  try:return all(isinstance(x,str) and datetime.date.fromisoformat(x).isoformat()==x for x in v)
  except (ValueError,TypeError):return False
 if kind=='clearance':
  return (exact(v,['mode','value_m','sigma_m']) and v['mode']=='measured' and number(v['value_m']) and number(v['sigma_m'])) or (exact(v,['mode','verified']) and v['mode']=='on_bed_verified' and v['verified'] is True)
 if kind=='discharge_sigma':return exact(v,['sigma_m3_s','definition']) and number(v['sigma_m3_s']) and text(v['definition'])
 if kind=='ice_bounds':return exact(v,['lower_C','upper_C']) and number(v['lower_C'],minimum=-273.15,maximum=0) and number(v['upper_C'],minimum=v['lower_C'],maximum=0)
 if kind=='limits':return exact(v,['mass_relative','energy_relative']) and number(v['mass_relative']) and number(v['energy_relative'])
 if kind=='grid':
  if not exact(v,['spacings_m','event_times_s','absolute_error_s','method']):return False
  dx=v['spacings_m'];tt=v['event_times_s']
  return isinstance(dx,list) and isinstance(tt,list) and len(dx)>=3 and len(dx)==len(tt) and all(number(x,minimum=1e-300) for x in dx) and all(a>b for a,b in zip(dx,dx[1:])) and all(number(x) for x in tt) and number(v['absolute_error_s']) and text(v['method'])
 return False

def valid_evidence(e):
 if not exact(e,['source_uri','sha256','locator','derivation','status']):return False
 uri=e['source_uri']
 try:
  u=urlparse(uri) if isinstance(uri,str) else None
  good_uri=bool(u and u.scheme=='https' and u.hostname and not u.username and not u.password and not any(c.isspace() for c in uri))
  if good_uri:u.port  # invalid port syntax must refuse, not crash
  if good_uri and '[' in u.netloc:ipaddress.IPv6Address(u.hostname)
 except (ValueError,TypeError,UnicodeError):return False
 return good_uri and sha(e['sha256']) and text(e['locator']) and text(e['derivation']) and e['status'] in ('observed','derived','verified_numerical')

def check(meta):
 issues=[]
 if not isinstance(meta,dict):return dict(schema_version=2,metadata_complete=False,may_issue_environmental_prediction=False,issues=['root must be an object'],scientific_review_required=True)
 if type(meta.get('schema_version')) is not int or meta.get('schema_version')!=2:issues.append('schema_version must be integer2')
 for g in sorted(set(meta)-set(RULES)-{'schema_version'},key=str):issues.append('unknown group '+str(g))
 values={}
 for group,rules in RULES.items():
  data=meta.get(group)
  if not isinstance(data,dict):issues.append(group+': missing object');continue
  for key in sorted(set(data)-set(rules),key=str):issues.append(group+'.'+str(key)+': unknown field')
  for key,kind in rules.items():
   name=group+'.'+key;record=data.get(key)
   if not exact(record,['value','evidence']):issues.append(name+': typed value/evidence record required');continue
   if not valid_value(kind,record['value']):issues.append(name+': invalid '+kind)
   else:values[name]=record['value']
   if not valid_evidence(record['evidence']):issues.append(name+': auditable evidence record required')
 limits=values.get('numerics.acceptance_limits_predeclared')
 if limits:
  for field,bound in [('mass_relative_residual','mass_relative'),('energy_relative_residual','energy_relative')]:
   value=values.get('numerics.'+field)
   if value is not None and value>limits[bound]:issues.append('numerics.'+field+': exceeds supplied predeclared limit')
 return dict(schema_version=2,metadata_complete=not issues,may_issue_environmental_prediction=False,scientific_review_required=True,issues=issues,
  interpretation='Metadata completeness only. Source authenticity, model applicability, calibration independence and predictive accuracy require scientific review; they are not certified here.')

def template():
 return dict(schema_version=2,**{g:{k:dict(value=None,evidence=dict(source_uri=None,sha256=None,locator=None,derivation=None,status=None)) for k in fields} for g,fields in RULES.items()})

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('metadata',type=Path,nargs='?');p.add_argument('--template',type=Path);a=p.parse_args()
 if a.template:
  a.template.write_text(json.dumps(template(),indent=2)+'\n');return
 if not a.metadata:p.error('supply metadata path or --template PATH')
 try:
  def reject_constant(x):raise ValueError('nonfinite JSON literal '+x)
  def unique(pairs):
   d={}
   for k,v in pairs:
    if k in d:raise ValueError('duplicate JSON key '+k)
    d[k]=v
   return d
  raw=a.metadata.read_bytes();m=json.loads(raw,parse_constant=reject_constant,object_pairs_hook=unique);result=check(m);result['input_sha256']=hashlib.sha256(raw).hexdigest()
 except (ValueError,TypeError,OverflowError,RecursionError,OSError) as exc:result=dict(metadata_complete=False,may_issue_environmental_prediction=False,issues=[str(exc)],scientific_review_required=True)
 print(json.dumps(result,indent=2,allow_nan=False));raise SystemExit(0 if result['metadata_complete'] else 2)
if __name__=='__main__':main()
