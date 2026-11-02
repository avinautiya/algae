import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from validity_contract_v2 import RULES,check,template

def documented_fixture():
 """Deliberately synthetic; complete syntax must NEVER certify field validity."""
 d=template();e=dict(source_uri='https://example.org/synthetic-fixture',sha256=hashlib.sha256(b'SYNTHETIC NOT FIELD EVIDENCE').hexdigest(),locator='synthetic test record',derivation='only tests schema parsing',status='derived')
 vals={'dates':['2026-01-01'],'nonnegative':.01,'text':'documented synthetic convention','true':True,'clearance':dict(mode='measured',value_m=.1,sigma_m=.01),'fraction':.9,'discharge_sigma':dict(sigma_m3_s=.01,definition='source standard uncertainty'),'contact':'observed','ice_bounds':dict(lower_C=0,upper_C=0),'sha256':e['sha256'],'grid':dict(spacings_m=[.1,.05,.025],event_times_s=[10.,11.,11.2],absolute_error_s=.5,method='synthetic grid study'),'flags':[True,True],'limits':dict(mass_relative=.1,energy_relative=.1)}
 for g,fields in RULES.items():
  for k,kind in fields.items():d[g][k]=dict(value=copy.deepcopy(vals[kind]),evidence=copy.deepcopy(e))
 return d

class ContractTests(unittest.TestCase):
 def test_template_refuses(self):self.assertFalse(check(template())['metadata_complete'])
 def test_complete_never_claims_valid_science(self):
  r=check(documented_fixture());self.assertTrue(r['metadata_complete']);self.assertFalse(r['may_issue_environmental_prediction']);self.assertTrue(r['scientific_review_required'])
 def test_string_booleans_refuse(self):
  for v in ['false','true',0,1,[],{},None,False]:
   d=documented_fixture();d['geometry']['width_measured_at_sensor_section']['value']=v;self.assertFalse(check(d)['metadata_complete'])
 def test_uncertainty_not_bool_negative_or_nan(self):
  for v in [True,-1,float('nan'),float('inf'),'0.1']:
   d=documented_fixture();d['geometry']['floor_survey_sigma_m']['value']=v;self.assertFalse(check(d)['metadata_complete'])
 def test_fraction_domain(self):
  for v in [-.1,1.1,True]:
   d=documented_fixture();d['temperature']['observed_support_fraction']['value']=v;self.assertFalse(check(d)['metadata_complete'])
 def test_numerics_required(self):
  d=documented_fixture();del d['numerics'];self.assertFalse(check(d)['metadata_complete'])
 def test_unsettled_hydraulics(self):
  d=documented_fixture();d['numerics']['hydraulic_convergence_per_interval']['value']=[True,False];self.assertFalse(check(d)['metadata_complete'])
 def test_residual_exceeds_limit(self):
  d=documented_fixture();d['numerics']['mass_relative_residual']['value']=.2;self.assertFalse(check(d)['metadata_complete'])
 def test_evidence_is_required(self):
  for value in [{},[],None,dict(source_uri='unknown',sha256='0'*64,locator='none',derivation='unknown',status='assumed')]:
   d=documented_fixture();d['measurement']['independent_of_calibration_data']['evidence']=value;self.assertFalse(check(d)['metadata_complete'])
 def test_schema_and_unknown_typos(self):
  for change in ['schema','field','group']:
   d=documented_fixture()
   if change=='schema':d['schema_version']=True
   elif change=='field':d['geometry']['widht']=True
   else:d['numeric']={}
   self.assertFalse(check(d)['metadata_complete'])
 def test_clearance_and_ice_bounds(self):
  d=documented_fixture();d['geometry']['sensor_clearance_m_or_on_bed_verified']['value']=True;self.assertFalse(check(d)['metadata_complete'])
  d=documented_fixture();d['contact']['ice_temperature_known_or_bounded']['value']=dict(lower_C=-1,upper_C=1);self.assertFalse(check(d)['metadata_complete'])
 def test_cli_duplicate_and_nonfinite_rejection(self):
  with tempfile.TemporaryDirectory() as tmp:
   path=Path(tmp)/'bad.json'
   for content in ['{"schema_version":2,"schema_version":2}','{"x":NaN}']:
    path.write_text(content);p=subprocess.run([sys.executable,str(Path(__file__).with_name('validity_contract_v2.py')),str(path)],capture_output=True,text=True);self.assertEqual(p.returncode,2);self.assertFalse(json.loads(p.stdout)['metadata_complete'])
if __name__=='__main__':unittest.main()
