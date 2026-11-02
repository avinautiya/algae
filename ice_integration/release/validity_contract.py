"""Version2 adapter for the release metadata gate.

`check` keeps the legacy result keys, but never grants scientific authority.
Exit0 means typed metadata completeness ONLY; exit2 means refused/incomplete.
Legacy inputs require explicit migration; there is no permissive fallback.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
if __package__:
 from . import validity_contract_v2 as core
else:
 import validity_contract_v2 as core

MAX_INPUT_BYTES=8*1024*1024

def load_strict(path):
 if Path(path).stat().st_size>MAX_INPUT_BYTES:raise ValueError('metadata exceeds8MiB limit')
 raw=Path(path).read_bytes()
 if len(raw)>MAX_INPUT_BYTES:raise ValueError('metadata exceeds8MiB limit')
 def constant(x):raise ValueError('nonfinite JSON literal '+x)
 def unique(pairs):
  result={}
  for key,value in pairs:
   if key in result:raise ValueError('duplicate JSON key '+key)
   result[key]=value
  return result
 return json.loads(raw,parse_constant=constant,object_pairs_hook=unique),hashlib.sha256(raw).hexdigest()

def check(meta):
 result=core.check(meta)
 unmet={}
 for issue in result['issues']:
  head=issue.split(':',1)[0]
  group,key=head.split('.',1) if '.' in head and head.split('.',1)[0] in core.RULES else ('_schema',head)
  if key not in unmet.setdefault(group,[]):unmet[group].append(key)
 result['unmet']=unmet
 result['label_if_not']='conditional component hindcast/scenario; metadata incomplete, not field validated' if not result['metadata_complete'] else 'metadata complete; scientific review required, not field validated'
 return result

def migrate_legacy(meta,input_sha256=None):
 """Retain every legacy value; do not turn assertions into evidence.

 Wrongly typed values remain wrong and must be corrected from source, never
 converted. Unknown fields are preserved in the separate migration record.
 """
 if not isinstance(meta,dict):raise ValueError('legacy root must be an object')
 if 'schema_version' in meta:raise ValueError('explicit migration accepts unversioned legacy input only')
 result=core.template();mapped=[];unknown=[]
 for group,fields in meta.items():
  if group not in core.RULES or not isinstance(fields,dict):
   unknown.append(dict(group=group,value=copy.deepcopy(fields)));continue
  for key,value in fields.items():
   if key not in core.RULES[group]:unknown.append(dict(group=group,key=key,value=copy.deepcopy(value)));continue
   result[group][key]['value']=copy.deepcopy(value)
   mapped.append(dict(field=group+'.'+key,value_type_matches=core.valid_value(core.RULES[group][key],value),evidence_imported=False))
 report=dict(migration_version=2,input_sha256=input_sha256,legacy_document=copy.deepcopy(meta),mapped=mapped,unknown_fields=unknown,
  interpretation='Values copied without coercion. Evidence remains blank; all missing numerical records remain blank. Nothing here verifies a source or makes a false/unknown field true.')
 return result,report

def write_new(path,value):
 path=Path(path)
 if path.exists():raise ValueError('refusing to overwrite existing output '+str(path))
 path.parent.mkdir(parents=True,exist_ok=True)
 with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,delete=False) as f:
  name=Path(f.name)
  try:json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
  except Exception:name.unlink(missing_ok=True);raise
 try:
  # link creates destination exclusively; an intervening writer is not replaced.
  os.link(name,path)
 finally:name.unlink(missing_ok=True)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('metadata',type=Path,nargs='?')
 p.add_argument('--template',type=Path)
 p.add_argument('--migrate',type=Path,help='explicitly migrate an unversioned legacy document')
 p.add_argument('--output',type=Path)
 p.add_argument('--migration-report',type=Path)
 a=p.parse_args()
 if sum(x is not None for x in [a.metadata,a.template,a.migrate])!=1:p.error('choose metadata, --template, or --migrate')
 try:
  if a.template:
   if a.output or a.migration_report:p.error('--template takes no other outputs')
   write_new(a.template,core.template());return
  if a.migrate:
   if not a.output or not a.migration_report:p.error('--migrate requires --output and --migration-report')
   paths=[x.resolve() for x in [a.migrate,a.output,a.migration_report]]
   if len(set(paths))!=3:raise ValueError('source and both outputs must be distinct')
   if a.output.exists() or a.migration_report.exists():raise ValueError('migration outputs already exist; refusing overwrite')
   source,digest=load_strict(a.migrate);migrated,report=migrate_legacy(source,digest)
   # Both files deliberately remain incomplete evidence; never issue approval.
   write_new(a.migration_report,report);write_new(a.output,migrated)
   print(json.dumps(dict(status='migrated_unverified',metadata_path=str(a.output),report_path=str(a.migration_report),may_issue_environmental_prediction=False)))
   return
  if a.output or a.migration_report:p.error('output options require --migrate')
  value,digest=load_strict(a.metadata);result=check(value);result['input_sha256']=digest
 except (ValueError,TypeError,OverflowError,RecursionError,OSError) as exc:
  result=dict(schema_version=2,metadata_complete=False,may_issue_environmental_prediction=False,scientific_review_required=True,issues=[str(exc)],unmet={'_schema':['invalid input']},label_if_not='refused input; not field validated')
 print(json.dumps(result,indent=2,allow_nan=False))
 raise SystemExit(0 if result['metadata_complete'] else 2)

if __name__=='__main__':main()
