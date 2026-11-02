# Independent final result audit

**Early-time analytic accuracy is not resolved:** FV64 with400steps per diffusion time has11.78–16.41% cumulative ice-heat error at0.1tau. Its<0.055%error at1tau is a different regime and does not certify0.1tau. All12cases show much larger change under temporal200→400 refinement at32cells than spatial32→64 refinement at400steps. Do not call the early FV reference converged.

## Recomputed checks

```json
{
  "coupled_records": 30,
  "all30complete": true,
  "frozen_source_hashes_match": true,
  "frozen_cases_hash_match": true,
  "frozen_spec_hash_match": true,
  "all_input_hashes_match": true,
  "all_record_freeze_hashes_match": true,
  "all_configurations_match": true,
  "contrast_count": 18,
  "maximum_recomputed_contrast_difference": 0.0,
  "all_fixed_footprint_gates_pass": true,
  "max_abs_phase_geometry_volume_residual_m3": 4.3021142204224816e-15,
  "max_record_relative_energy_residual": 2.2541140484893632e-15,
  "boundary_records": 12,
  "maximum_recomputed_boundary_difference": 0.0,
  "boundary_early_error_range_percent": [
    11.783667740166216,
    16.407917887099586
  ],
  "boundary_1tau_max_abs_error_percent": 0.0540140149503987,
  "boundary_10tau_max_abs_error_percent": 0.0001885880708929477,
  "maximum_geometry_summary_total_difference": 0.0,
  "boundary_time_refinement_dominates_space_at_early_time_all12": true
}
```

## Findings

- **FR1 — major numerical qualification:** At0.1diffusion time FV64/400steps is NOT an accurate analytic reference:11.78–16.41% cumulative ice-heat error. At fixed32cells, doubling steps200→400 reduces error much more than32→64cells at400steps in all12cases. Early-time accuracy requires explicit temporal refinement; late1tau error<0.055% cannot be extended to early time.
- **FR2 — scope:** Coupled reference is same governing implementation with refined settings, not independent CFD or field truth. Contrasts establish approximation differences only.
- **FR3 — scope:** Covered-footprint area is matched in frozen gates. Total-face and floor geometries evolve; geometric-study floor melt is not an identical-Eulerian-footprint accuracy comparison. Separate whole-model output differences from matched-cover contrasts.
- **FR4 — scope:** Prescribed equal wet/dry net heating is a controlled boundary intervention, not matched natural weather; noexchange control is not a professional glacier baseline. Geometric undercut has adiabatic ledges and no collapse/lifetime inference.
- **FR5 — figure interpretation:** Figure01 draws small11–16% FV early error on same yaxis as333–656% RC error; data are correct but FV error can look visually negligible. Label early FV residual numerically or add separate inset before presenting accurate early-time reference. Figure02uses correctly computed matchedcovered differences; Figure03directly plots saved heightbands and honestly labels adiabatic ledges/no collapse.

All30coupled and12boundary outputs are retained. All contrast arithmetic, analytical error arithmetic, configuration/source hashes, fixed-covered-footprint gates and geometry-summary totals were checked independently without editing parent source, freezes, results or figures. No new Nature-level novelty or environmental superiority is inferred.

Source inspected:run_study.py,run_boundary_study.py,analyze_study.py,study/SPEC.md,study/FREEZE.json,study/COUPLED_CASES.json and all savedcase JSON. Figure calculations were traced to source; no figure rewrite or expensive rerun was performed.
