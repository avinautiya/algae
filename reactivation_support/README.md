# Independent support for the reactivation milestone

This package speeds verification and analysis without changing Claude's active solver or runs. Source audit is pinned to commit c70c55319232d1dae3d3a7ee856327101df995ea; later Claude edits may already supersede a finding.

Start with **DIAGNOSTIC_FINDINGS.md**. The most useful lead is a matched inlet-geometry test that reduces the initial partition discrepancy substantially. The package also provides an ODE-free prismatic reference (~416× faster in the measured reference benchmark), phase-exchange verification, exact event scheduling, connectivity diagnostics and campaign analysis.

## Use from repository root

```sh
python3 -m unittest discover -s reactivation_support -t . -p 'test_*.py'
python3 -m reactivation_support.verify_snapshot
python3 -m reactivation_support.geometry_audit
python3 -m reactivation_support.partition_diagnostic --output diagnostics/partition.json
python3 -m reactivation_support.smooth_control_diagnostic
```

Core helper tests: 46. These verify the helper implementation and its assumptions; they do not pass the existing solver's full verification gate. Diagnostic snapshot runs are development evidence only.

For completed campaign files (replace the example input path):

```sh
python3 -m reactivation_support.campaign_analysis results/strip.json --threshold 0.001 --persistence-hours 24 --output diagnostics/strip_summary.json
python3 -m reactivation_support.plot_campaign results/strip.json --output figures/strip --label 'Development run — not final evaluation'
```

Declare the threshold and persistence interval in the frozen evaluation specification rather than tune them on results. If flow is pulsed, supply `--ordinary-from-day` explicitly to prevent flood-only wetting from being called ordinary-flow recovery.

## Integration priorities

1. Read the matched entrance diagnostics before buying more refinement for the unexplained partition discrepancy. Use an applicable reference; preserve and quantify the abrupt-step case instead of changing production geometry to force agreement.
2. Replace the expensive GVF ODE reference only where the prismatic assumptions apply. Do not reuse it as an exact 2-D strip truth.
3. Use exact event-aligned morphology intervals and log actual removal times, convergence flags and storage.
4. Choose a consistent phase-mass/energy formulation. The local exchange adapter is opt-in; it does not solve the multirate cumulative ledger by itself.
5. Use the analysis and plotting tools once the independently verified, frozen campaign finishes. Core analysis does not need a full simulation rerun.

## Evidence and scope

- `snapshot/`: read-only text snapshot of the three pinned source files used for reproduction, not a production fork.
- `records/`: independently generated diagnostics and reference benchmarks.
- `records/example_development_run.json`: short deliberately non-final run exercising output compatibility and removal-time quantization.
- `figures/`: diagnostic/development figures only.
- `phase_update.py`: conservative exchange with an explicit sensible-energy contract.
- `event_schedule.py`: event-aligned morphology intervals.
- `hydraulic_reference.py`: exact integrated flat-channel GVF length relation, with scalar root inversion.
- `campaign_analysis.py`: interval/event/uncertainty diagnostics without fabricated final-data labels.
- `connectivity.py`: geometric minimax path, not a hydraulic discharge model.

The source pin, normalization and local file hashes are recorded in SNAPSHOT_PROVENANCE.json. No historical originality, new glacier discovery or improved environmental forecast is asserted for these debugging tools.
