"""Reproduce the declared CN calculation and expose unresolved target authority.

No density or calendar convention is chosen by its agreement with the paper.
This module does not convert SR50 lowering into validated pure-ice melt.
"""
import argparse,hashlib,importlib.util,json,math
from pathlib import Path

def error_percent(model_mwe,lowering_m,rho_ice):
 if not all(type(x) in (int,float) and math.isfinite(x) and x>0 for x in [model_mwe,lowering_m,rho_ice]):raise ValueError('positive finite quantities required')
 return 100*(model_mwe/(lowering_m*rho_ice/1000)-1)

def density_sensitivity(d):
 model=d['model_sum_m_water_equivalent'];lowering=d['thresholded_surface_lowering_m'];published=d['published_Table3_percent']
 raw_target=model/(1+published/100)
 rows=[dict(rho_kg_m3=rho,error_percent=error_percent(model,lowering,rho)) for rho in [900.,915.,917.]]
 return dict(protocol_density_kg_m3=900.,protocol_density_authority='original FINAL_RUNS protocol',
  protocol_density_role='prescribed conversion of MODEL OUTPUTS; it does not independently document observed-SR50 postprocessing',
  sensitivity_only=rows,
  rho915_authority='final supplement MCC19 individual-model definition; NOT an observation-processing rule',
  target_implied_by_published_rounded_error_mwe=raw_target,
  implied_target_is_not_an_observation=True,
  density_that_exactly_fits_rounded_score_kg_m3=1000*raw_target/lowering,
  fitting_density_to_score_is_forbidden=True,
  published_score_reproduced_under_authorized_conditional_rule=d['published_score_reproduced'],
  independent_melt_validation=False,
  unresolved=['exact final processed observation total','daily interval alignment/aggregation','ice-versus-snow and surface-displacement interpretation','calibration versus independent evaluation provenance'])

def execute(original_code,protocol_dir,data_zip,outputs):
 spec=importlib.util.spec_from_file_location('pinned_reproduction',original_code);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 d=module.calculate(protocol_dir,data_zip,outputs)
 return dict(conditional_reproduction=d,target_authority_audit=density_sensitivity(d),
  executed_reproduction_source_sha256=hashlib.sha256(original_code.read_bytes()).hexdigest(),
  result='NO_GO for independent pure-ice validation; no alternative convention selected by fit')

def main():
 p=argparse.ArgumentParser();p.add_argument('--original-code',type=Path,required=True);p.add_argument('--protocol-dir',type=Path,required=True);p.add_argument('--data-zip',type=Path,required=True);p.add_argument('--outputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 if a.out.exists():raise ValueError('output exists')
 result=execute(a.original_code,a.protocol_dir,a.data_zip,a.outputs);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2));raise SystemExit(2)
if __name__=='__main__':main()
