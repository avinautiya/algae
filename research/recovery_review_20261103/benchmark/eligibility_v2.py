"""Claim-specific benchmark gates; metadata cannot certify scientific accuracy."""
import json

COMMON=['model_output_member_crc_verified','model_output_sha256_recorded','actual_thickness_column_verified','model_version_mapping_verified','date_axis_verified','matching_observation_window_verified','units_and_sign_verified','snow_target_definition_verified','paper_metric_reproduced']
INDEPENDENT=['parameter_calibration_excludes_test_targets','candidate_frozen_before_scoring','support_frozen_before_scoring','issue_time_feature_provenance_verified','evaluation_not_selected_for_model_success']

def check(meta,*,claim):
 if claim not in ['published_baseline_reproduction','independent_pure_ice_melt_test']:raise ValueError('unsupported claim')
 if not isinstance(meta,dict):raise ValueError('metadata object required')
 required=COMMON+(INDEPENDENT+['pure_ice_melt_observation_verified'] if claim=='independent_pure_ice_melt_test' else [])
 unmet=[k for k in required if meta.get(k) is not True]
 # This is the paper's final-cumulative-melt comparison; daily MAE alone does
 # not reproduce its Table3 ranking or percentage error.
 if meta.get('metric')!='final_cumulative_percentage_error':unmet.append('metric must match paper final cumulative percentage error')
 return dict(claim=claim,metadata_gate_passed=not unmet,unmet=unmet,accuracy_improvement_established=False,interpretation='A passed metadata gate allows a qualified comparison; it does not authenticate assertions or establish independence, field accuracy, novelty or superiority.')

if __name__=='__main__':
 import argparse
 from pathlib import Path
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('metadata',type=Path);p.add_argument('--claim',choices=['published_baseline_reproduction','independent_pure_ice_melt_test'],required=True);a=p.parse_args();r=check(json.loads(a.metadata.read_text()),claim=a.claim);print(json.dumps(r,indent=2));raise SystemExit(0 if r['metadata_gate_passed'] else 2)
