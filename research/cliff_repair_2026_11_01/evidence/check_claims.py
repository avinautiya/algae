"""Recompute corrected claims; fail closed for the unavailable full-v4 ROI gate.

This checks a pinned evidence package, not scientific truth or novelty. Legacy
mechanism numbers are report-derived and authenticated only to their Git blobs,
not reclassified here from raw traces or reproduced by a new solver campaign.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def git_blob(data):
 return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()

def evaluate(registry,root=ROOT):
 issues=[];checked=[]
 if not isinstance(registry,dict):return dict(evidence_consistent=False,issues=['registry must be an object'],may_claim_full_v4_matched_melt=False,may_claim_field_superiority=False)
 manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text())
 for entry in manifest['files']:
  p=root/'source'/entry['path']
  if not p.is_file() or git_blob(p.read_bytes())!=entry['blob_sha']:issues.append('source mismatch: '+entry['path'])
  else:checked.append(entry['path'])
 if issues:return dict(evidence_consistent=False,issues=issues,may_claim_full_v4_matched_melt=False,may_claim_field_superiority=False)
 if type(registry.get('schema_version')) is not int or registry.get('schema_version')!=1 or registry.get('audited_commit')!=manifest['ref']:issues.append('registry version/ref mismatch')
 facts=registry.get('facts',{})
 if not isinstance(facts,dict):issues.append('facts must be an object');facts={}
 numbers=json.loads((root/'source/ice_integration/release/figures/ERROR_BUDGET_NUMBERS.json').read_text())
 errors=list(numbers['test_signed_errors_nominal_m'].values());bias=sum(errors)/len(errors);mae=sum(abs(x) for x in errors)/len(errors)
 weights=numbers['calibration_influence']['weight'];leverage=1/sum(x*x for x in weights)
 expected=dict(nominal_test_count=len(errors),nominal_pooled_bias_m=bias,nominal_pooled_MAE_m=mae,linear_imputer_marginal_bias_over_27_scenarios_m=numbers['factor_effects']['imp']['linear']['mean_pooled_bias_m'],leverage_concentration_equivalent_count=leverage)
 for key,value in expected.items():
  actual=facts.get(key)
  if type(actual) not in (int,float) or not math.isfinite(actual) or not math.isclose(actual,value,rel_tol=1e-12,abs_tol=1e-12):issues.append('incorrect recomputed fact: '+key)
 if len(errors)!=6 or not all(x>0 for x in errors) or bias!=mae:issues.append('nominal signed-error basis changed')
 v=facts.get('mechanism_versions',{})
 if not isinstance(v,dict):issues.append('mechanism_versions must be an object');v={}
 # These expected values are explicitly report-derived, not fresh simulations.
 expected_versions=dict(current_bed_operator='v4',raw_rescore_removal_classes=dict(maintained=3,genuine_reconnection=1,failed_by_horizon=2),distinct_pre_loss_histories=2,solid_strip_loss_comparison_histories=1,contact_persistence_boundary_days=[10.0,10.5],strength_defined_boundary_days_approx=9.2,day10_removal_reconnects_then_loses_contact_at_day=12.4,legacy_locked_ROI_melt_relative_errors=[-.54,-.67],legacy_locked_ROI_melt_absolute_errors_m3=[-.002,-.0034],v3_v4_matched_scorer_refusals_each=10,complete_v4_matched_ROI_campaign_done=False)
 for key,value in expected_versions.items():
  if v.get(key)!=value:issues.append('incorrect qualified report fact: '+key)
 report=(root/'source/ice_cliff_survival/reactivation/REPORT_M3_RESCORE.md').read_text()
 for token in ['3 maintained / 1 genuine reconnection / 2 failed','10.0 and 10.5','9.2 d','12.4 d','one solid-strip history','−54 %','−67 %','a v4 re-run of the full evaluation grid']:
  if token not in report:issues.append('pinned report support changed: '+token)
 for key in ['loss_GCI_is_statistical_CI','finest_grid_shift_is_not_full_coupled_error_bound']:
  desired=key.startswith('finest')
  if v.get(key) is not desired:issues.append('uncertainty meaning changed: '+key)
 if facts.get('leverage_count_is_not_statistical_independence_ESS') is not True:issues.append('leverage is not statistical ESS')
 gates=registry.get('gates',{})
 if not isinstance(gates,dict):issues.append('gates must be an object');gates={}
 for key in ['full_v4_matched_melt_ready','field_superiority_established','historical_originality_established','independent_field_validation']:
  if gates.get(key) is not False:issues.append('unsupported claim enabled: '+key)
 claims=registry.get('claims',[])
 if not isinstance(claims,list) or not all(isinstance(c,dict) for c in claims):issues.append('claims must be an object list');claims=[]
 ids=[c.get('id') for c in claims]
 if ids!=[f'C{i:02}' for i in range(1,18)]:issues.append('claim inventory must retain all17claims')
 statuses={c.get('id'):c.get('status') for c in claims if isinstance(c.get('id'),str)}
 for key,value in dict(C05='false',C08='unavailable',C13='qualified',C16='unsupported',C17='unestablished').items():
  if statuses.get(key)!=value:issues.append('claim ceiling violated: '+key)
 return dict(evidence_consistent=not issues,issues=issues,source_blobs_verified=len(checked),recomputed=expected,
  may_claim_full_v4_matched_melt=False,may_claim_field_superiority=False,historical_originality_established=False,
  v4_gate_reason='No authenticated v4 matched-ROI scorer/complete campaign supplied. Legacy ROI scores cannot satisfy this gate.',
  limitations='Report numbers are checked against pinned source text; this does not authenticate legacy solver runs, establish a joint causal error budget, or certify originality.')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--registry',type=Path,default=ROOT/'CLAIM_REGISTRY.json');p.add_argument('--target',choices=['v4-matched-melt','evidence-consistency'],default='v4-matched-melt');p.add_argument('--out',type=Path);a=p.parse_args()
 try:result=evaluate(json.loads(a.registry.read_text()))
 except (ValueError,KeyError,TypeError,OSError) as e:result=dict(evidence_consistent=False,issues=[str(e)],may_claim_full_v4_matched_melt=False,may_claim_field_superiority=False)
 if a.out:a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps(result,indent=2,allow_nan=False))
 passed=result['evidence_consistent'] and (a.target=='evidence-consistency' or result['may_claim_full_v4_matched_melt'])
 raise SystemExit(0 if passed else 2)
if __name__=='__main__':main()
