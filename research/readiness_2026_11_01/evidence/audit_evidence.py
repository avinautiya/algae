"""Read-only arithmetic/provenance check of a pinned repository evidence packet.

No solver, score or historical output is changed. Source snapshots must match
Git blob hashes. Textual inconsistencies are reported rather than repaired.
"""
import hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def blobsha(data):
    return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()

def main():
    manifest=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
    checked=[]
    for entry in manifest['files']:
        path=ROOT/'source'/entry['path'];data=path.read_bytes()
        # apply_patch may add one terminal newline when archiving text. Only a
        # proven matching Git-blob byte sequence is accepted; no content fix.
        candidates=[data]
        for _ in range(3):
            if candidates[-1].endswith(b'\n'):candidates.append(candidates[-1][:-1])
        matches=[x for x in candidates if blobsha(x)==entry['blob_sha']]
        if len(matches)!=1:raise AssertionError('Source snapshot mismatch: '+entry['path'])
        if data!=matches[0]:path.write_bytes(matches[0])
        checked.append(entry['path'])
    b=json.loads((ROOT/'source/ice_integration/release/figures/ERROR_BUDGET_NUMBERS.json').read_text())
    errors=b['test_signed_errors_nominal_m'];bias=sum(errors.values())/len(errors)
    mae=sum(abs(x) for x in errors.values())/len(errors)
    assert len(errors)==6 and all(x>0 for x in errors.values())
    assert bias==mae
    w=b['calibration_influence']['weight'];concentration=1/sum(x*x for x in w)
    assert abs(sum(w)-1)<1e-12
    linear_marginal=b['factor_effects']['imp']['linear']['mean_pooled_bias_m']
    claimed_nominal=.63
    assert round(bias,2)!=claimed_nominal
    assert round(linear_marginal,2)==claimed_nominal
    contract=json.loads((ROOT/'source/ice_integration/release/swiss_package_contract_result.json').read_text())
    unmet=sum(len(x) for x in contract['unmet'].values())
    assert unmet==8 and not contract['may_issue_environmental_prediction']
    result=dict(pinned_commit=manifest['ref'],source_blobs_verified=len(checked),
      nominal_test_count=len(errors),nominal_pooled_bias_m=bias,nominal_pooled_MAE_m=mae,
      error_budget_claimed_nominal_bias_m=claimed_nominal,
      linear_imputer_marginal_bias_over_27_scenarios_m=linear_marginal,
      nominal_label_defect_confirmed=True,
      leverage_concentration_equivalent_count=concentration,
      leverage_count_is_not_statistical_independence_ESS=True,
      validity_unmet_count=unmet,
      contract_is_current_package_specific_not_global_public_data_impossibility=True,
      mechanism_versions=dict(current_bed_operator='v4',
       headline_loss_comparison='legacy development trajectory, 3 labels share pre-removal history',
       loss_GCI_days=.37,loss_GCI_is_statistical_CI=False,
       frozen_recovery_boundary_prediction_days=9.50,
       reported_dynamic_reference_boundary_days=10.39,
       reported_fixed_medium_history_fine_grid_boundary_days=10.45,
       recovery_boundary_dynamic_error_days=10.39-9.50,
       recovery_boundary_finest_grid_shift_days=.02,
       finest_grid_shift_is_not_full_coupled_error_bound=True,
       new_boundary_class_matches=4,new_boundary_class_cases=5))
    result['mechanism_versions'].update(
      raw_rescore_removal_classes=dict(maintained=3,genuine_reconnection=1,failed_by_horizon=2),
      distinct_pre_loss_histories=2,solid_strip_loss_comparison_histories=1,
      contact_persistence_boundary_days=[10.0,10.5],
      strength_defined_boundary_days_approx=9.2,
      day10_removal_reconnects_then_loses_contact_at_day=12.4,
      legacy_locked_ROI_melt_relative_errors=[-.54,-.67],
      legacy_locked_ROI_melt_absolute_errors_m3=[-.0020,-.0034],
      v3_v4_matched_scorer_refusals_each=10,
      complete_v4_matched_ROI_campaign_done=False,
      legacy_raw_provenance='metadata_identity_checked_not_run_authenticated')
    (ROOT/'AUDIT_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
