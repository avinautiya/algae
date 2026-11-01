# Additive repair handoff — preserve the frozen M3 campaign

Import this folder separately; do not overwrite your active source or existing
results. Run `python3 verify.py`. It runs41 tests and launches no CFD campaign.

## Immediate source check before spending finest-grid compute

Inspect the current `Sim.lower_bed` against pinned `snapshot/swe2d.py:175`.
The pinned global `maximum(zb-dz,pool_floor)` mutates non-ice tail geometry even
at dt=0 and can count latent melt without equivalent geometric lowering.
Run the independent reproducer and review `bed_audit/FINDINGS.md`. If already
fixed on your branch, document the fixing commit and rerun the preflight there;
do not reapply a competing change blindly.

The supplied `bed_audit/operator.py` implements ice-only actual lowering.
An optional ice-domain floor causes explicit refusal rather than hidden clipping.
Removing the old floor is a versioned model change and requires inlet/outlet,
storage and phase accounting checks. It is not a complete liquid-mass/momentum
or multirate energy repair. A genuine bedrock contact requires its own thermal
closure. Do not simply clip the bed and keep consuming latent heat.

The paired coarse diagnostic preserves early loss with a.0237d time shift, but
this is not permission to transfer the old.37d GCI to the revised model. Re-run
the relevant grid sequence after the corrected operator is verified/refrozen.

## Repair the existing melt comparison without rerunning CFD first

Run from this folder:

    python3 -m melt.repair_eval /ACTUAL/results/eval snapshot/frozen/production_spec.json --output /SEPARATE/matched_melt.json

This uses existing initial/final `toe_floor_x` profiles over x3.5..4.5m across
the full0.5m toe width. It compares exactly the same physical ROI with supported
reduced controls, retaining original frozen predictions. Fifteen tests cover
scope, profile precision, identity/horizon checks and safe output paths.

Check every full-reach profile-versus-cumulative disagreement; do not assume
heat-based melt equals actual ice lowering where the old floor clipped.
The5e-6m³ central-ROI profile precision bound is neither model error nor a CI.
Save corrected errors and the whole/control/open-reach decomposition. Do not
assert improved melt accuracy until those actual scores support it.

The raw campaign JSONs are absent from the public pinned tree, so export them
and a run/config/source manifest. Raw results do not contain full Q/T provenance;
matching filenames or hashes alone cannot authenticate the execution setup.

## Correct recovery and independence claims

Use `classification/recovery.py` on full-resolution raw exports, supplying
actual removal/ordinary-return times. Keep theta and one-day persistence frozen.
Require sustained time AFTER assessment begins; a final wet sample is not proof.
Future exports should preserve both pre- and post-event diagnostics because
the current records only store the post-removal state at the event timestamp.

Replace the phrase 'four recovery outcomes' with separate maintained contact,
genuine reconnection and failed recovery counts. Report-derived nonporous
counts are2/1/1; all removal-case counts are3/1/2. Confirm against raw traces.
E04/E06/E07 share one pre-loss history. Do not use those identical loss times
as three independent validations, bootstrap units or three independent test cases.

Preserve the frozen scalar event times separately from the downsampled traces.
Those non-bare traces stop at day20; don't extrapolate to21d or substitute their
daily interpolated crossing for the original quarter-day event estimate.
If producing augmented traces, retain the frozen artifact and record a new
deterministic export version; evaluate the actual terminal model state rather
than append a repeated last sample.

## Refine only after those gates pass

`numerics/` reproduces the old conditional GCI and supplies a bounded isolated
runner. A3.5d fine run tests early loss only, not day12 removal. Aday12 recovery
test needs at least13d to observe one-day persistence; reproducing nine days of
absence requires the full21d. Use the same verified model and initial condition
on every grid, and explain the changing divider-thickness approximation.

Keep the old campaign as historical evidence, revised solver results as a new
version, and raw-data scoring repairs distinct from new model computations.
