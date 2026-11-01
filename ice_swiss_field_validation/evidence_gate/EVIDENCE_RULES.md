# Rules for the resumed Swiss comparison

These rules are written before the resumed scoring campaign, after the paper,
notebook outputs and some survey values have already been examined. This is a
retrospective evaluation specification, not a blind or externally preregistered
study.

1. Use original numerical survey and temperature records recovered from the
   author-disclosed ETH mirror. Preserve original member names and byte hashes.
   The mirror archive's identity with the later DOI deposit is unverified.
2. The primary physical target is change of directly surveyed channel-floor
   elevation, in metres, over exactly specified survey intervals. A waterline
   mark, fitted Nusselt number, or modelled melt volume is not that target.
3. Use the same accepted intervals, clock convention, geometry assumptions and
   supplied forcings for every model. Keep missing forcing as missing; do not
   bridge a long gap or fit missing inputs from the target incision.
4. Keep standard uncalibrated correlations, a training-only constant melt-rate
   predictor, and a training-only constant heat-transfer predictor. A calibrated
   physical model must be compared against equally data-informed baselines.
5. Calibrate only on the specified earlier intervals. Do not choose a diameter,
   clock shift, temperature offset, model family or exclusion from later errors.
   Report all frozen scenarios. Existing author-derived Nu from incision cannot
   be used as an independent feature for predicting that same incision.
6. Compare errors for each station as well as pooled errors. P3 and P5 share a
   channel, forcing and drainage event. They do not establish two independent
   environmental replications, much less superiority on Himalayan glaciers.
7. Keep measurement uncertainty, engineering sensitivity and statistical
   uncertainty separate. Source survey standard deviations are not automatic
   hard bounds. If +/-0.2 m is used as a conservative endpoint scenario, label
   that assumption; it is not a 95% confidence interval.
8. Use timestamp-weighted temperature integrals for the native-record analysis.
   Preserve author sample-weighted means as a separate version. Identify missing
   samples, calibration offsets, averaging windows and clock uncertainty.
9. Near-floor measured temperature is not demonstrated bulk-water temperature.
   A heat-transfer coefficient fitted against it is an effective conditional
   coefficient, not automatically a validated turbulent transport constant.
10. Observed future water temperature, stage or discharge produces a forced
    hindcast. It does not test an advance forecast of those forcings or validate
    the complete hydraulic/debris/cliff model.

`paired_evidence.py` implements a conditional worst-case absolute-error advantage
for supplied measurement and prediction intervals. It does not manufacture error
bars. Positive advantage in every tested group is a within-data consistency
check; negative groups must remain visible even if the pooled score improves.

The desired outcome is a physically useful predictor with lower held-out error.
Whether that outcome is supported is an experimental result, not an assumption.
