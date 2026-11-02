import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from validity_contract import check,migrate_legacy
from test_contract import documented_fixture

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('audited_legacy_contract',ROOT/'legacy/validity_contract.py')
legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)

class AdapterTests(unittest.TestCase):
 def test_original_seven_exploits_reproduced_then_rejected(self):
  base={g:{k:True for k in ks} for g,ks in legacy.REQUIRED.items()}
  edits=[None,('measurement','independent_of_calibration_data','false'),('temperature','observed_support_fraction',2.),('geometry','floor_survey_sigma_m',-.1),('discharge','discharge_uncertainty',float('nan')),('geometry','cross_section_width_source',{}),('geometry','floor_elevation_survey_dates',[])]
  for edit in edits:
   value=copy.deepcopy(base)
   if edit:value[edit[0]][edit[1]]=edit[2]
   self.assertTrue(legacy.check(value)['may_issue_environmental_prediction'])
   # Even a source-shaped V2 wrapper must reject the malformed typed values.
   migrated,_=migrate_legacy(value)
   r=check(migrated);self.assertFalse(r['metadata_complete']);self.assertFalse(r['may_issue_environmental_prediction'])
 def test_old_result_keys_preserved_but_authority_removed(self):
  r=check(documented_fixture());self.assertTrue(r['metadata_complete']);self.assertEqual(r['unmet'],{});self.assertFalse(r['may_issue_environmental_prediction']);self.assertIn('scientific review',r['label_if_not'])
 def test_swiss_migration_no_coercion_or_fabrication(self):
  value=json.loads((ROOT/'legacy/swiss_package_metadata_status.json').read_text());value['new_group']={'unknown':[]};value['geometry']['new_field']='unverified'
  m,report=migrate_legacy(value)
  for group,fields in value.items():
   if group=='new_group':continue
   for key,x in fields.items():
    if key=='new_field':continue
    self.assertEqual(m[group][key]['value'],x)
    self.assertTrue(all(v is None for v in m[group][key]['evidence'].values()))
  self.assertEqual(report['legacy_document'],value);self.assertEqual(len(report['unknown_fields']),2)
  self.assertIs(m['geometry']['floor_survey_sigma_m']['value'],True)
  self.assertIs(m['measurement']['independent_of_calibration_data']['value'],False)
  self.assertFalse(check(m)['metadata_complete'])
 def test_no_automatic_legacy_fallback(self):
  base={g:{k:True for k in ks} for g,ks in legacy.REQUIRED.items()}
  self.assertFalse(check(base)['metadata_complete'])
  with self.assertRaises(ValueError):migrate_legacy(documented_fixture())
 def test_huge_integer_and_malformed_urls_refuse_without_crash(self):
  d=documented_fixture();d['geometry']['floor_survey_sigma_m']['value']=10**400
  self.assertFalse(check(d)['metadata_complete'])
  for uri in ['https://[not-an-ip]/x','https://example.org:badport/x','https://user:secret@example.org/x','https://exa mple.org/x']:
   d=documented_fixture();d['geometry']['floor_survey_sigma_m']['evidence']['source_uri']=uri;self.assertFalse(check(d)['metadata_complete'])
 def test_nested_json_types_never_approve_or_crash(self):
  bad=[None,{},[],[{}],{'a':[]},'unknown',True]
  fixture=documented_fixture()
  for group,fields in fixture.items():
   if group=='schema_version':continue
   for key in fields:
    for b in bad:
     d=copy.deepcopy(fixture);d[group][key]=copy.deepcopy(b);self.assertFalse(check(d)['metadata_complete'])
 def cli(self,*args):return subprocess.run([sys.executable,str(ROOT/'validity_contract.py'),*map(str,args)],capture_output=True,text=True)
 def test_cli_nonfinite_duplicate_deep_json(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'x.json'
   for raw in ['{"schema_version":2,"schema_version":2}','{"x":NaN}','{"x":Infinity}','['*2000+']'*2000]:
    p.write_text(raw);r=self.cli(p);self.assertEqual(r.returncode,2);self.assertFalse(json.loads(r.stdout)['may_issue_environmental_prediction']);self.assertNotIn('Traceback',r.stderr)
 def test_cli_explicit_migration_and_overwrite_protection(self):
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);m=d/'v2.json';report=d/'migration.json';source=ROOT/'legacy/swiss_package_metadata_status.json'
   args=['--migrate',source,'--output',m,'--migration-report',report]
   r=self.cli(*args);self.assertEqual(r.returncode,0);self.assertEqual(json.loads(r.stdout)['status'],'migrated_unverified');self.assertFalse(check(json.loads(m.read_text()))['metadata_complete'])
   before=m.read_bytes();r=self.cli(*args);self.assertEqual(r.returncode,2);self.assertEqual(m.read_bytes(),before)
   r=self.cli('--migrate',source,'--output',source,'--migration-report',d/'other.json');self.assertEqual(r.returncode,2)
 def test_complete_cli_zero_means_only_metadata_complete(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'complete.json';p.write_text(json.dumps(documented_fixture()));r=self.cli(p);self.assertEqual(r.returncode,0);d=json.loads(r.stdout);self.assertTrue(d['metadata_complete']);self.assertFalse(d['may_issue_environmental_prediction'])
 def test_repository_package_module_entrypoint(self):
  import shutil
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);release=d/'ice_integration/release';release.mkdir(parents=True)
   (d/'ice_integration/__init__.py').write_text('');(release/'__init__.py').write_text('')
   for name in ['validity_contract.py','validity_contract_v2.py']:shutil.copyfile(ROOT/name,release/name)
   meta=d/'fixture.json';meta.write_text(json.dumps(documented_fixture()))
   r=subprocess.run([sys.executable,'-m','ice_integration.release.validity_contract',str(meta)],cwd=d,capture_output=True,text=True)
   self.assertEqual(r.returncode,0,r.stderr);self.assertFalse(json.loads(r.stdout)['may_issue_environmental_prediction'])
if __name__=='__main__':unittest.main()
