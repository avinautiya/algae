# Explicit gap reconstruction for a sensitivity version

The strict native-temperature adapter is unchanged. `gap_reconstructed_temperature.py` is a separate, explicitly reconstructed forcing component. Its primary linear bridge was selected using the controlled **temperature masking** audit, where it beat both endpoint-persistence alternatives across all six predefined station/duration groups. No incision or melt target selected the method. Both persistence alternatives remain mandatory sensitivities; this is not an independent new field dataset or a claim of melt-model superiority.

Methods and support rules are retained in `GAP_RECONSTRUCTION_SPEC.json`, frozen before any subsequent melt rescore. A gap is reconstructed only between successive **finite valid measured endpoint values** no more than 24 hours apart. The original consecutive-finite ≤60-second support remains measured support; finite endpoint gaps spanning NaNs or longer acquisition intervals become reconstructed support. Gaps exceeding 24 hours remain unsupported in full, and there is no extrapolation outside the first/last finite values. Reconstructing missing values assumes the endpoints remained representative across sensor handling or changing immersion; the masking audit cannot guarantee that assumption.

```python
from gap_reconstructed_temperature import GapReconstructedTemperatureSeries
series = GapReconstructedTemperatureSeries.from_archive(
    original_subset_root, 5,
    max_observed_gap_seconds=60,
    max_reconstruction_gap_seconds=86400)
result = series.integrate_positive(
    start_naive_datetime, end_naive_datetime,
    extra_offset_C=0,
    method='linear_endpoint_bridge')
```

Alternative methods are `constant_previous_endpoint` and `constant_next_endpoint`. General `.integrate()` is signed; `.integrate_positive()` exactly integrates `max(T+shift,0)` and retains the signed result separately. The extra offset is a declared sensitivity, not an additional source calibration. P3/P5 final +0.035/+0.040°C corrections are applied exactly once when loading original raw-suffix data.

Each call reports:

- Total, observed-support and reconstructed-support temperature integrals separately, in °C·s.
- Observed/reconstructed seconds and fractions of the requested window.
- Remaining unsupported fraction and longest unsupported interval.
- Reconstruction method, clamp flag, sensitivity shift and unresolved source clock.

The point-floor-to-bulk temperature offset remains unknown. The source clock remains naive; callers must declare any UTC/CEST alignment scenario. The 24-hour reconstruction limit is a methodological choice, not an accuracy guarantee. For gaps longer than six hours, the masking audit showed appreciable diurnal bias even for linear bridges.

Five added tests verify separate observed/reconstructed constant and linear integrals, persistence alternatives, no >24-hour fill or extrapolation, exactly-once source calibration and exact positive-temperature zero crossing. Together with the earlier tests, 15 tests pass.
