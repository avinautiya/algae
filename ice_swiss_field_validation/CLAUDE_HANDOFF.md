# Claude handoff — original Swiss data and completed comparisons

This package is additive. No active Claude worktree, checkpoint, frozen recovery campaign or solver source was modified. Do not replace existing solver results with these component hindcasts.

## Completed work and actual verdict

Original numerical data are now available; missing-archive access is resolved through the author-disclosed ETH mirror. Use `source/` and its manifest, not graph digitization as primary forcing. Retain the source version distinction from the later DOI deposit.

| Analysis | Actual outcome |
|---|---|
| Original native support, no reconstructed gaps | Every complete manual survey interval contains missing approved temperature support. No fully observed-forcing score. |
| Author-daily reconstruction | Scaled flow-sensitive DB improves fitted constant-h nominal MAE by 20.6% at P3 and 14.1% at P5, but loses to smooth Gnielinski at P5. |
| Native-time reconstruction, gaps <=24h | Same six model families; P5 failure persists across all three reconstructors and engineering scenarios. |
| Local stage-based geometry | On identical training/test cohorts, pooled MAE improves 23.6% versus its own global-geometry family; it still loses to simple baselines. Geometry availability changes training membership, so unmatched version scores are not a feature ablation. |
| Hydraulic sign correction | Source algebraic error confirmed twice. Corrected constant-friction velocity model fails its later-data comparison; correction is not proof of predictive superiority. |
| Temperature reconstruction | Linear bridging beats endpoint persistence in all six predeclared station-duration groups. This validates only a retrospective reconstruction component on masked observed blocks. |

**Do not claim consistent field superiority, new validated turbulent physics, full cliff melt-volume accuracy, or Himalayan transfer from these results.** All P3/P5 cases belong to one Swiss stream event. Models are conditional on imposed future observed/reconstructed forcing.

## Use the code concretely

- `temperature_digitization/native_temperature.py`: native time integration, exact positive-temperature clipping, source verification, and separate limited-gap reconstruction. Original strict coverage remains available. Calibration offsets apply exactly once; misleading `_corr_raw` filenames still need the final P3/P5 offsets.
- `gap_reconstruction_audit/`: frozen independent forcing-reconstruction test, all block predictions and figures. Do not turn its errors into guaranteed bounds on actual maintenance gaps.
- `conditional_benchmark/`: all V2–V5 freezes, source/code hashes, exact executed source snapshots, saved predictions, exclusions, independent audits, figures and `REPRODUCE.md`.
- `original_mirror/`: safe range recovery, direct source/head-expression evidence, nominal sensitivity and clock audit. Original author products are unmodified.
- `velocity_transfer/`: frozen transfer protocol, all seven observations, trained coefficients, four later test cases, and negative results.
- `evidence_gate/`: explicit paired-error consistency accounting; engineering ranges are not confidence intervals.

Run each folder's tests in an isolated copy and retain all result versions. Existing tests should pass without OpenFOAM. Plotting needs Matplotlib/NumPy; the numerical adapters and score accounting use the standard library. The independent arithmetic audits already pass to roundoff.

## Required integration behavior

1. Add a field-forcing adapter to the coupled model using the supplied temporal integration. Keep valid observed support, reconstructed support and unsupported time separate in outputs/restarts.
2. Do not infer zero flow/contact from a purged sensor record. Pressure head equals full channel depth only when bed clearance is known; see `OBSERVATION_OPERATOR.md`.
3. Distinguish true supplied geometry from pressure-proxy/reconstructed geometry. Recompute area, hydraulic diameter and Reynolds number consistently at each supplied section. A P3 Reynolds series must not silently stand in for P5 local hydraulics.
4. Allocate the water's finite incoming sensible heat only once across floor and face. A water-temperature-driven component hindcast is not a prediction of that future water temperature or discharge.
5. Carry the negative comparisons into the paper. A 23.6% improvement within one family is not superiority over the best baseline. Do not select a diameter, clock shift, temperature correction or imputer by melt test error.
6. Keep source and model uncertainty separate; a floor sensor is not an independent bulk-water thermometer, and source standard deviations cannot be silently relabelled as hard bounds or confidence intervals.
7. Any new model developed after these scores needs a new evaluation claim. Preserve these cases as inspected development evidence, not newly blind validation.

An independent Rochemelon2005 external-test feasibility task is being prepared separately. Its forcing, geometry, target and model-freeze gates must pass before any external accuracy claim. Do not substitute modelled/fitted source curves for observations or reuse Rochemelon-fitted coefficients as independent predictions on that event.
