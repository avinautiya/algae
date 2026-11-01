# Frozen external thermal prediction adapter

The adapter is executable. It contains **no Rochemelon bed observations, fitted site coefficient, built-in site geometry or field score**. Current matching does not establish the unlabeled Figure6 target's exact surveyed section or a supported thermal path from lake temperature to that point. That unresolved input gate must not be replaced with an invented channel length/depth.

## Run and verify

Python3.9+ standard library only. Keep the protocol/freezes beside the script.

```sh
python -m unittest discover -s model -p 'test_*.py'
python model/predict_external.py model/SYNTHETIC_PACKET.json --out synthetic_prediction.json
```

The13 checks cover analytic exponential cooling, identical-reach splitting, finite heat bounds, zero inlet heat, cold-ice sensible/latent accounting, bed/wall/volume partition, a stable tiny-NTU limit, verified critical point controls, data/clock refusal, uncapped statistical-reference identity, omitted-gravity diagnostic behaviour, and derived numerical geometry refusal. They are verification, not field validation.

`SYNTHETIC_PACKET.json` is explicitly synthetic. Its dimensions and fluxes must not be copied into a Rochemelon input packet. The executable stamps it synthetic and refuses promotion of synthetic forcing/geometry to an external field packet.

## Required packet contract

- `purpose`: `synthetic_verification` or `conditional_external_hindcast`.
- `clock_convention`: `UTC` with explicit timestamp offsets, or consistently `declared_source_naive` with a separately audited source-clock scenario.
- `alignment`: true declared `path_alignment_verified` and `geometry_independent_of_bed_target`, plus a nonempty source citation. These are evidence declarations; the parent source audit must actually establish them.
- `records`: contiguous positive-duration forcing intervals. Each has `start`, `end`, positive `Q_m3_s`, nonnegative `inlet_temperature_C`, explicit nonpositive `ice_temperature_C`, `forcing_source`, and `forcing_basis` (`interval_observed` or explicitly `sourced_interval_scenario` for external packets).
- Each record has the complete ordered thermal-route `sections`. Every section requires a unique `name`, positive `width_m`, explicit `length_m`, geometry `source`, and `geometry_basis` (`measured` or `source_assumption`). Uniform reaches need independently sourced positive `depth_m` and `length_m`.
- A `point_control` has `length_m:0` and either explicit sourced local depth or `depth_mode:critical_control`. The latter also needs `critical_control_verified:true` and `critical_control_source`; its depth follows critical-flow physics from Q/width. It is a local flux target, not a finite-volume reach. Unknown upstream thermal distance or depth is still a refusal.
- Optional sourced `potential_drop_m` with `potential_drop_source` is diagnostic only; it supplies no assumed friction heat.

CLI exits2 with `status:REFUSED` when the numerical/provenance contract is incomplete. No temporal gap filling, extrapolation, default field width/depth, geometry inference from observed melt, or model fitting is performed. A legitimate source-assumption scenario remains conditional rather than measured. Actual source quality is not proven by checking a JSON flag.

## What the output means

All seven frozen model outputs are retained. Six have finite-water-heat physical component ledgers; the original uncapped constant-rate statistical reference stays separate and warns if its requested heat exceeds supply. H-based models cool water exponentially through uniform imposed reaches; the heat-limited constant-rate physical variant cools linearly until its supply is exhausted. Downstream points receive the cooled water, not unchanged lake temperature.

The output separates watts/joules, latent heat, fresh-ice sensible cost, bed/sidewall energy and volume, mean finite-footprint bed lowering, and zero-area point bed lowering. Imposed hydraulic Q is not updated with produced meltwater; the produced mass fraction is explicitly reported. Atmosphere and dissipation are omitted for all methods. Cold-ice heat cost is a declared algebraic sensitivity, **not** a cold-ice conduction PDE. Tube-correlation numerical domains do not establish rough open-channel validity.

After the parent input/target matching gate passes, prepare a source-backed forcing packet and run the same script to obtain component predictions. Scoring is a separate operation on fixed matched observations. No result from this module alone establishes independent field accuracy or consistent superiority.
