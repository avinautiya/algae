import copy,json,subprocess,sys,tempfile,shutil,unittest
from pathlib import Path
from check_claims import evaluate,ROOT
from build_docs import documents

class ClaimsTests(unittest.TestCase):
 def registry(self):return json.loads((ROOT/'CLAIM_REGISTRY.json').read_text())
 def test_corrected_consistency_without_promoting_claims(self):
  r=evaluate(self.registry());self.assertTrue(r['evidence_consistent'],r);self.assertEqual(r['source_blobs_verified'],11);self.assertFalse(r['may_claim_full_v4_matched_melt']);self.assertFalse(r['may_claim_field_superiority'])
 def test_nominal_not_marginal(self):
  r=self.registry();r['facts']['nominal_pooled_bias_m']=.63;self.assertFalse(evaluate(r)['evidence_consistent'])
 def test_shared_history_not_three_validations(self):
  r=self.registry();r['facts']['mechanism_versions']['solid_strip_loss_comparison_histories']=3;self.assertFalse(evaluate(r)['evidence_consistent'])
 def test_uncertainty_not_ci(self):
  r=self.registry();r['facts']['mechanism_versions']['loss_GCI_is_statistical_CI']=True;self.assertFalse(evaluate(r)['evidence_consistent'])
 def test_contact_strength_endpoints_not_collapsed(self):
  r=self.registry();r['facts']['mechanism_versions']['strength_defined_boundary_days_approx']=10.4;self.assertFalse(evaluate(r)['evidence_consistent'])
 def test_missing_v4_data_cannot_enable_gate(self):
  r=self.registry();r['gates']['full_v4_matched_melt_ready']=True;self.assertFalse(evaluate(r)['evidence_consistent']);self.assertFalse(evaluate(r)['may_claim_full_v4_matched_melt'])
 def test_unsupported_superiority_and_novelty_stay_blocked(self):
  for key in ['field_superiority_established','historical_originality_established','independent_field_validation']:
   r=self.registry();r['gates'][key]=True;self.assertFalse(evaluate(r)['evidence_consistent'])
 def test_claim_status_cannot_be_promoted(self):
  r=self.registry();r['claims'][12]['status']='verified';self.assertFalse(evaluate(r)['evidence_consistent'])
 def test_source_tamper_refuses(self):
  with tempfile.TemporaryDirectory() as temp:
   p=Path(temp);shutil.copy2(ROOT/'SOURCE_MANIFEST.json',p/'SOURCE_MANIFEST.json');shutil.copytree(ROOT/'source',p/'source');q=p/'source/ice_integration/release/ERROR_BUDGET.md';q.write_text(q.read_text()+'\nchanged\n');self.assertFalse(evaluate(self.registry(),root=p)['evidence_consistent'])
 def test_generated_documents_are_exact(self):
  for name,content in documents(self.registry()).items():self.assertEqual((ROOT/'proposed_release'/name).read_text(),content)
 def test_default_cli_fails_closed_but_consistency_can_pass(self):
  script=str(ROOT/'check_claims.py');p=subprocess.run([sys.executable,script],capture_output=True,text=True);self.assertEqual(p.returncode,2);self.assertFalse(json.loads(p.stdout)['may_claim_full_v4_matched_melt']);p=subprocess.run([sys.executable,script,'--target','evidence-consistency'],capture_output=True,text=True);self.assertEqual(p.returncode,0)
 def test_malformed_nested_types_refuse(self):
  for bad in [None,[],True]:self.assertFalse(evaluate(bad)['evidence_consistent'])
  for key in ['facts','gates','claims']:
   r=self.registry();r[key]=True;self.assertFalse(evaluate(r)['evidence_consistent'])
if __name__=='__main__':unittest.main()
