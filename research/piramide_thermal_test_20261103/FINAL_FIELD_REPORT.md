# Locked field campaign — actual findings

**The candidate does not outperform every declared baseline.** The first field period gave three nominal baseline wins; the subsequent unchanged-model batch does not beat the air-temperature index and is effectively tied with the buried-temperature index. Neither more solver tolerance nor choosing a favorable capacity/geometry scenario changes that conclusion.

## Three later periods, frozen before their outcomes

Training remains December2–16,2022. The candidate, derivative comparator, geometry/capacity scenarios, forcing reconstructions and all39 gains remain unchanged. Later observed ice ablation is59.046117,61.451154 and89.054069 mm water equivalent in the three predeclared consecutive periods throughFebruary10,2023.

| Predeclared method | Three later-period MAE (mm w.e.) | Four-period MAE including earlier result (mm w.e.) |
|---|---:|---:|
| Forward thermal candidate |14.220367|11.754572|
| Derivative-OLS thermal comparator |14.349710|12.378427|
| Constant rate |18.173057|20.480931|
| Buried-temperature index |14.192518|12.248945|
| Air-temperature index |9.870827|7.906105|

The later-period primary metric was declared as mean absolute period error before the later camera outcomes were opened. Total-error cancellation is not substituted for it. The four-period summary is secondary and retains the earlier air-index loss. The fourth proposed later period was excluded for every model before its observations were opened, because the air record endsFebruary21 and cannot support the February24 endpoint without extrapolation. Its camera outcomes remain unopened.

All39 fixed predictions were published at `ebd47166c8bd1715eff3668337d5f18d4ead4e52` before the later camera endpoints were released. This is a locked retrospective evaluation on published field data. The source paper was consulted, so it is not a prospective field trial or a guarantee that no participant knew any source-paper result.

## What was implemented and checked

The package now includes a source-verified MATLAB MCOS decoder, native temperature/moisture/weather tables, exact camera-unit binding, shared-reading uncertainty arithmetic, source-bounded geometry and capacity scenarios, conservative debris/ice conduction with a complementary melting interface, deterministic freezes, role-separated observation access, prediction seals, fixed-gain transfer and independent recomputation.

Version1's negative-basal-heat physical no-go is preserved. Version2's missing endpoint is preserved as unscorable. Version3's nearest-finite replacement was selected solely by source time/availability, with every model and gain held fixed and revised predictions published before releasing its value. That source-coverage amendment is explicitly disclosed.

The native camera displacement difference times917kg/m³ reproduces every scored target. Independently reading the supplied ablation column confirms that its endpoint differences are metres water equivalent, agreeing after multiplication by1000. The frozen operator was not changed to improve a score. Quantitative camera precision and internal correction/reset uncertainty remain unresolved.

Long-trajectory mesh/deep-boundary checks change primary period predictions by at most about0.0192% in the tested refinements. Energy residuals are at roundoff scale, complementarity holds and no negative latent melt is generated. Those checks support numerical correctness under assumptions; they do not explain or erase the field prediction errors.

## Scientific limits that affect the conclusion

This is one tower on one glacier, with serial periods sharing observation endpoints. It is an observed-temperature-conditioned hindcast, not a future-weather forecast. Calendar periods, scenarios and camera frames are not independent glacier replications. Camera uncertainty is not a source-certified confidence interval; conditional bounds are labelled as assumptions. All finite scenarios remain visible, rather than selecting the best after scoring.

The physical solver's energy ledger applies to its unscaled raw melt component. A training-only empirical output gain corrects that component; the scaled predictions are not independently measured heat capacity or newly verified physical latent-energy totals. In the primary case, raw training melt is about362.52mm w.e. against85.20mm observed, requiring a correction near0.235. That discrepancy is important evidence, not a detail to hide behind calibration.

Static geometry, initial ice temperature, neglected ice thermal advection and wet-debris phase/transport remain assumptions. Full professional energy-balance models are not reproduced here. The derivative comparator is OLS, not the source NTLS method. The results therefore cannot establish general superiority over professional glacier models or a novel glacier discovery.

## Further hypothesis gate

An additive contact-resistance study investigates whether imperfect debris–ice contact could explain the amplitude discrepancy. It is not a claimed discovery: porous-ice contact resistance and thermal-resistance networks already exist in the literature. A preliminary training-only quasi-steady check also argues against a single constant contact resistance explaining the entire factor-of-four discrepancy while preserving the measured near-base temperature. This remains a restricted falsifier with transient/storage qualifications, not a universal impossibility theorem.

The second tower's public data are eligible for a future independent-stake comparison, but its camera outcome values remain unopened. Any revised method needs a justified identification mechanism and a separately frozen test. The four scored periods are now evaluation evidence and cannot be reused as fresh validation of variants tuned to beat their errors.
