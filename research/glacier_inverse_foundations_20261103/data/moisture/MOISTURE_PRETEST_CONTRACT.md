# Moisture transfer audit before Tower2 camera release

Tower2 camera outcomes remain unopened. This report uses process observations and published weather inputs only. `MOISTURE_TRANSFER_FEATURE_AUDIT.json` and `audit_moisture_transfer.py` reproduce all figures below. Common half-open calendar periods use the Tower2exact03:00anchors; the previous Tower1Dec30endpoint was00:00, so these are feature summaries, not recomputed Tower1comparison errors.

## What the source measures

The [source paper](https://doi.org/10.1029/2025JF008360) places Tower1moisture sensors at0.30/0.40/0.50m under nominal0.50m debris; Tower2sensors at0.19/0.29/0.39m under0.45m debris. Each depth/thickness has±0.02m stated uncertainty. Thus SM3's nominal distance above ice is0m atTower1and0.06m atTower2; separate geometry boxes permit0–0.04m and0.02–0.10m respectively. They are not interchangeable interface-water measurements.

[Manufacturer specifications](https://www.onsetcomp.com/products/sensors/s-smc-m005) describe typical accuracy±0.031m³/m³, not3%of the reading, for0–50°C mineral-soil operation. [The manual](https://www.onsetcomp.com/sites/default/files/resources-documents/15081-J%20S-SMx%20Manual.pdf) warns of temperature dependence and operation outside the accurate range below freezing. [Calibration guidance](https://www.onsetcomp.com/sites/default/files/resources-documents/15921-C%20Calibration%20Equations%20for%20S-SMC-M005.pdf) requires material-appropriate conversion. This is a manufacturer engineering reference, not a field calibration certificate or guaranteed error bound for these rock-debris installations. No local soil-specific calibration was found in the inspected release.

The sensor samples a small surrounding volume. Near an ice boundary, changes in material/contact geometry can affect its dielectric response. Treating SM3as a direct ice-interface saturation or whole-column water fraction is unsupported. The reported native negative/out-of-porosity readings are retained and flagged, not silently repaired.

## Feature trajectories

| Period start | T1SM3mean | T2SM3mean |
|---|---:|---:|
| Dec2training | .050193 | .070419 |
| Dec16 | .031822 | .074917 |
| Dec30 | .036304 | .093503 |
| Jan13 | .025626 | .087390 |
| Jan27 | .017369 | .088346 |

Tower1dries in this proxy while Tower2becomes wetter. That contrast is useful for a genuinely unscored transfer test; it is not evidence that the same fitted relation will improve Tower2melt. All672lateTower1SM3readings are below0.031; the late mean is less than the manufacturer's typical absolute accuracy. Eleven Tower1SM3values are negative in the Dec30–Jan13period. Tower2SM3is finite and within0–assumedphi=.44 for all five periods.

The source weather member was selectively retrieved (77.5KB compressed) and ZIP-CRCverified. It contains336hourly rows per period, with supplied precipitation totals0/6.6000/24.0000/0/0mm. This is an ERA5-labelled calibration input with measured/filled/lapse-corrected provenance, not a directly observed local rain record. **All snow values and outgoing-shortwave values are missing** in this member. No direct local albedo or native snow-free flag is available here. The source model masks snow through a separate albedo product, which was not downloaded. Missing snow cannot become zero by assumption.

AtTower1the DT3/DT4counts below−0.4°C are0/0during training,65/2in Dec16,70/18in Dec30,17/5in Jan13, and0/0in Jan27. DT5never drops below−0.4°C. AtTower2DT3/4/5never drop below−0.4°C in these windows. This is consistent with warmer lower debris atTower2, but does not prove phase state at every moisture probe; no exactly co-located temperature exists atTower2SM3depth0.39m. Upper debris still freezes, so a full-column liquid-only model remains conditional.

## Normalizing by own training moisture

If a sensor obeys y=a*theta+b, then y/y_train cancelsa only whenb=0. Additive offsets distort relative factors, especially near zero. A stableb cancels in y−y_train, while multiplicative differences remain. Therefore relative normalization does not erase cross-site calibration or different sensing depth.

For illustration only, using±0.031as a separately bounded engineering error reference, the observed Tower1training-to-late decrease(.032824) is smaller than the possible two-reading error span(.062). Independent temporal errors therefore do not establish its direction. Under a stable additive offset the direction survives but its relative magnitude does not. Do not present either calculation as a confidence interval or validated sensor-error model.

A non-singular anomaly transform, e.g. exp(gamma*(y−mean_train_y)), can be evaluated as a **sensor-proxy predictor**. It is not automatically a calibrated physical k(theta)law. An intrinsic-conductivity interpretation additionally requires sensor calibration, spatial support, pore phase and mass/enthalpy consistency. Choosing a different sensor after seeing Tower2errors is prohibited.

## Pretest input and comparison contract

1. Freeze source/member/CSV hashes, Tower2SM3as the nominated proxy, its actual depth and geometry uncertainty, the common source-clock interpretation, all five intervals and the exact half-hour integration operator. Preserve the declared 18h maximum air-temperature gap rule.
2. Freeze the proxy transform and all hyperparameters using **Tower1development evidence only**. State that Tower1melt outcomes were repeatedly inspected during development. Do not call Tower1a fresh validation or tune gamma/probe/window against Tower2evaluation outcomes.
3. Permit one Tower2training-ablation gain per applicable predictor. One training target cannot identify two new free melt coefficients. Any additional moisture coefficient must be fixed using Tower1development data or independent physical/thermal fitting before Tower2camera release.
4. Give the air-temperature-plus-moisture baseline the **same extra moisture parameter budget, Tower1development access and tuning opportunity** as the gradient-plus-moisture candidate. Retain air-only, gradient-only and constant-rate comparators. Report every frozen model and all four periods, not the best scenario after testing.
5. Declare whether raw negative moisture values are statistical sensor-proxy inputs or physically invalid observations. A physical-water model must refuse invalid values or use a predeclared measurement-state treatment; clipping to0/phi is not an observation. Do not silently drop flux-time contributions. Tower2nominatedSM3currently passes its physical-value feature gate, but that is not a calibration certificate.
6. Retain the existing signed heat/cold-storage operator. A positive gradient or empirical dose is not automatically gross ice melt. If a reduced estimator uses an imposed latent-only boundary or a positive-part dose, label that as a conditional empirical operator, rather than claiming a conservative complete glacier model.
7. Fix the phase assumptions: sensor temperatures near0are uncertain, upper debris can freeze, pore-ice latent storage is not represented by finite sensibleC, and local snow state is missing. Do not infer snow absence from zero supplied precipitation or a wetness trend. State the model's supported regime before scoring.
8. Publish model code, selected parameters and normalized predictions before training-target release. Publish all issued final predictions before evaluation-target release. Bind source camera operator/density917and intervals into the release receipt. Source camera uncertainty and interior correction/reset provenance remain unresolved.
9. Score all four intervals with identical observation operators. Treat the four windows as repeated observations at one site, not four independent glaciers. A moisture term that beats a weak gradient comparator but loses to the equally informed air-plus-moisture model does not establish a superior predictor.

Moisture-dependent debris heat transfer is established prior art, including [Giese et al.2020](https://doi.org/10.5194/tc-14-1555-2020). The defensible new result would be a predeclared successful transfer or an interpretable failure under these gates; no historical novelty or performance outcome is presumed.
