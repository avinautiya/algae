# Tower2 contact-boundary feasibility gate

No Tower2 camera distance or ablation value has been exposed. Camera data was consumed only through native timestamp arrays and a finiteness mask. `TOWER2_FEATURE_ELIGIBILITY.json` records the selection policy, all timestamp/finiteness evidence, source SHA and feature coverage. `audit_tower2_features_only.py` reproduces it.

## Eligible untouched target windows

All six fixed nominal endpoints are exact finite native camera labels:2December03:00,16December03:00,30December03:00,13January03:00,27January03:00,10February03:00. No nearest-time substitution was needed. They define one training and four held-out intervals. All labels lie on the complete native1800-second thermal grid; no collision occurs. Native timezones are empty, so preserve the common local-source clock without claiming UTC.

The four held-out intervals have maximum native2m air-temperature gaps of16.167,15.333,2.333 and1.667hours; training maximum is3hours. All satisfy the predeclared18h gap ceiling and source bounds cover them. The long first two gaps remain substantial interpolation assumptions, not measured weather. Freeze the exact time-integration operator and sensitivity policy before camera release. Do not fill gaps using melt observations.

The final [source paper](https://doi.org/10.1029/2025JF008360), table1, reports Tower2's45cm debris, thermistors at0/11/22/33/44cm and moisture probes19/29/39cm; thickness/depth uncertainty±2cm. Tower2 falls12February, after the proposed final endpoint. Its nominal deepest sensor is1cm above the interface; location uncertainty is comparable to or larger than that gap. The nominal surface sensor is gravel-covered. These are not exact boundary locations.

## Training-only capacity and temperature evidence

For672half-hour samples in[2December03:00,16December03:00), all three moisture series satisfy0<=theta<=assumedphi=.44; none was clipped or discarded. Means are0.265067113,0.088673065 and0.070419494. Their unweighted mean0.141386558 is not a measured whole-column volume average.

Independent source lithology gives rock density2657.313333kg/m³ and specific heat816.126667J/kg/K. Using assumedphi=.44 and the fixed-phase mixture formula from the preceding audit gives conditional bulk sensible capacities:

- Dry:1.215005MJ/m³/K.
- Training probe proxy:1.805830MJ/m³/K.
- Saturated liquid pores:3.053674MJ/m³/K.

These are physical sensitivity scenes, not local measured capacity bounds or confidence intervals. Wet phase changes require enthalpy rather than arbitrarily enlarged sensibleC.

TrainingDT5 ranges0.034–1.310°C;393of672samples exceed0.4°C. Thus the nominal last probe is not directly demonstrated to be ice at the melting point. Sensor error, unobserved gap, cold-ice storage and contact resistance must not be collapsed into a claimed measuredR.

## What temperature data can identify

For constant k,C and a fixed zero-temperature ice boundary behind a contact resistanceR, the temperature equations depend on

alpha=k/C, and ell=kR.

Indeed k→a*k, C→a*C, R→R/a leaves both alpha andell invariant for any positivea. All temperature predictions remain unchanged, while conductive heat and inferred melt scale bya. Therefore temperature fitting alone cannot separately identify k,C andR under this model. An independently justifiedC or a separate heat/melt observation is required to fix the scale. This elementary structural scaling argument is not a new historical theorem.

The effective boundary length also absorbs an unsensed basal gapg: in a quasi-steady gap with the same conductivity, ell_effective=g+kR_contact. The nominalg=.01m is uncertain; separate±.02m thickness/depth boxes permit physicalg from0 to.05m. Do not report fittedell/k as intrinsic contact resistance without resolving that geometry.

An exploratory training-only temperature diagnostic fits T5≈ell*(T4−T5)/.11, yieldingell≈.076635m. This is **not** an independently establishedR estimate: predictor and response shareDT5 measurement error, and unsteady storage violates the quasi-steady approximation. It is a motivation for a predeclared dynamic-model comparison, not proof of a contact mechanism or its originality.

## Fair next test

Fit temperature-supported alpha/effective boundary length using training temperatures only, carry the declared geometry/C scenes, and compare a nonnegative-contact model with its R=0 restriction and the already strong air-temperature baseline. Use the same single training-ablation gain allowance for every applicable model, if such calibration is declared. Freeze every model and gain rule before releasing training endpoints, publish final predictions before all four evaluation observations, and report all four windows including failures.

Tower2 is a new target location with unopened camera values, but shares glacier/weather and source-processing methods with Tower1. It tests transfer of an idea developed after Tower1's failures; it does not create independent glaciers or establish Nature-level novelty. No camera precision or manual-reset error bound is presently justified.
