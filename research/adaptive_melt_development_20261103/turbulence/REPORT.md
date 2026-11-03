# Meteorological forcing and neutral-gradient turbulent-flux audit

**New usable source inputs are available.** The small original per-site forcing files supply pressure and radiation-related meteorology, and the observed four-level profiles permit a conditional neutral-gradient H/LE feature calculation. These are not independently validated turbulent fluxes or a complete surface energy balance. No new melt target was read and no model parameter was fitted to melt.

## Source and geometry gate

The original [Melo-Velasco et al.2025 paper](https://doi.org/10.1029/2025JF008360), section2.2.1/Table1, and the archive Readme agree on heights0.5,1,2,2.7m for levels1–4. T/RH accuracies are0.35C/3%; wind accuracies0.45m/s at the lower levels and0.89m/s above. Towers1/2 remain upright through the selected Dec2–Feb10 window. The paper describes pressure derived from DGA stations and ERA5-Land adjustment rather than a local tower barometer. The source profile method already estimates aerodynamic parameters; this is not new turbulence physics.

The source [Zenodo archive](https://doi.org/10.5281/zenodo.11581073) includes per-site `run_data_Piramide_Tower1/2...mat` inputs of279KB each. Selected ranged downloads were size/CRC32 verified and SHA256 recorded; whole-archive MD5 verification is not claimed. `export_verified_forcing.py` uses an explicit meteorological-variable whitelist, not camera or optimized-parameter arrays. Each site has4188 hourly records. `Pre` is hPa, confirmed by the source `Aerodynamic_Resistence.m` conversion toPa. `Latm` feeds the incoming-longwave variable in the source driver; `ea`/`esat` arePa. Raw solar partitions are retained individually. No value for pressure or instrument height was invented.

### Clock correction to the initial concern

The MAT starts at Nov1T00 while the CSV starts at Nov1T01, but this is a coverage difference, not evidence of a1h time shift. At4187 matching naive labels, Tower1 MATPre minus CSV pressure is constant−32.415939326hPa with standard deviation1.02e−8hPa; shifting labels by±1h raises that deviation to0.342hPa. More directly, the original TC driver converts MATLABdatenum and selects `dateTimes==AWS.DATE`, while the calibration script parses the CSV labels and retimes source temperatures by their labels. Therefore equal naive labels reproduce the author's binding, with no target-dependent lag selection. The model input has `DeltaGMT=-4`; source tower timezone fields remain empty, so externally verified physical UTC alignment is still not claimed. `t_bef=1.5,t_aft=-0.5` are preserved solar-integration metadata, not silently replaced by a symmetric window.

## Implemented neutral-profile equations

At each native common tower label, fit ordinary least squares against ln(z) at the four verified heights:

    U(z)=a_u+b_u ln(z), T(z)=a_T+b_T ln(z), q(z)=a_q+b_q ln(z).

Take u*=kappa*b_u, kappa0.4. For positive b_u, neutral constant-flux, zero-displacement, equal scalar diffusivity assumptions give upward-positive

    H_up = -rho*c_p*kappa*u* * b_T
    LE_up = -rho*L_v*kappa*u* * b_q.

Here q is specific humidity computed from the source RH/T and actual site pressure using q=.622e/(P−.378e); rho=P/[R_d*T*(1+.61q)]. Constants are explicit in `neutral_profile.py`. Surface temperature and saturation at the debris surface are not required for these gradient proxies. This does not eliminate uncertainty in atmospheric stability, flux divergence, roughness-sublayer structure, displacement height or sensor offsets.

The approximate gradient Richardson diagnostic at1m is g*b_Tv/(meanTv*b_u²). It is reported, not used to silently delete unstable records. Neutrality is not established by a fitted logarithmic profile alone. H_up and LE_up must be subtracted in a downward-net surface energy budget. They must not also be subtracted from an observed conductive gradient as though they were an additional independent loss; that would double-count atmospheric effects already expressed in the measured thermal profile.

### Explicit refusals and uncertainty

A nonpositive wind log slope produces an unsupported-row status and empty flux fields, never zero heat. The wind-slope uncertainty diagnostic is the sum of absolute OLS coefficients times the independently reported per-level accuracies; the analogous T-slope radius uses0.35C. These are conditional sensitivity diagnostics, not Gaussian confidence intervals. No accuracy is divided by sqrt(number of samples). RH calibration, height uncertainty and stability are not exhausted by these two bounds.

## Actual target-free coverage

| Diagnostic | Tower1 | Tower2 |
|---|---:|---:|
| Native profile rows in Dec2–Feb10 window |9469|8838|
| Nominal supported neutral proxies |9465|8834|
| Explicit nonpositive-wind-slope refusals |4|4|
| Wind-slope signs resolved under stated accuracy |2210|2854|
| Temperature-slope signs resolved under stated accuracy |3669|4781|
| Rows with nominal absolute gradientRi<=0.03 |3798|4325|
| Largest native gap |16.83h|16.17h|

Nominal H ranges are−90.1..344.7W/m² and−157.0..405.3W/m²; nominal LE ranges−152.7..173.4W/m² and−143.1..287.7W/m². These are diagnostic ranges, not measured or bounded actual flux. The large missingness and weakly resolved gradients mean full-period native turbulent-energy coverage does not pass a strict gate. A separate explicitly reconstructed scenario can be defined before comparison, but cannot be called a fully observed surface-energy budget.

The four refusal timestamps and native-gap inventory are in `FEATURE_COVERAGE.json`. Retain these flags when a downstream code builds period exposures; dropping unsupported rows or turning blanks into zero is forbidden.

## Artifacts and tests

- `Tower1_MODEL_FORCING.csv`, `Tower2_MODEL_FORCING.csv`: source pressure, temperature, radiation components and other whitelisted meteorological inputs, including raw MATLABdatenum.
- `Tower1_NEUTRAL_PROFILE_FEATURES.csv`, `Tower2_NEUTRAL_PROFILE_FEATURES.csv`: native labelled H/LE proxies with source and model-validity diagnostics.
- `FORCING_EXPORT_MANIFEST.json`, `NEUTRAL_FEATURE_MANIFEST.json`, `ALL_MEMBERS_MANIFEST.json`: hashes, units, raw member provenance and explicit scope.
- `neutral_profile.py`, `build_neutral_features.py`: executable feature pipeline.
- `test_neutral_profile.py`: five tests pass for exact logarithmic profiles, zero scalar flux, upward warm/moist transport signs, inversion refusal and physical input validation.

Scripts were read, not executed as MATLAB. No camera archive, optimized conductivity, optimized roughness or source model melt output was loaded. Radiation/pressure inputs are independently processed meteorology, not independent new local measurements.

## Remaining gate for a physical predictor

This unlocks a substantially better forcing basis than air temperature alone. It does not by itself resolve radiometric surface temperature (DT1 is buried), net longwave/absorbed shortwave at the pit, ice cold content, latent water transport or the source's vertical flux-constancy assumptions. A state-space closure can use these inputs with explicit storage and boundaries, but field superiority must be tested. The existing failed predictions remain unchanged. No Nature significance or historical originality is inferred from finding these source files.

## Final source-code binding addendum

The CRC-verified `MAIN_FRAME.m`, lines434–436, explicitly defines total solar input as STOT=SAB1+SAB2+SAD1+SAD2 before spectral/directional albedo weighting. `ShortwaveFluxesOS.m` documents incoming shortwave inW/m² and separately computes absorbed/reflected components. The exported partition sum is therefore incoming shortwave, not already absorbed. `Surface_Temperature_Debris.m` uses the surface residual Rnet−H−QE−G+Qv. These source statements support the signs and operator binding; they do not establish a measured pit albedo, surface temperature or validity of the neutral-gradient closure.

`ALL_MEMBERS_MANIFEST_V2.json` is the complete current eight-member size/CRC32/SHA manifest. Earlier `FETCH_MANIFEST.json` and `ALL_MEMBERS_MANIFEST.json` are retained as historical intermediate records. `ARCHIVE_INDEX.json` is copied into this packet and bound to SHA256149537e63f200115c655de7aaa779483ade2e252d5c1587cca9272adb8ad34f1. The fetcher now uses this local metadata index, verifies its hash and writes a new exclusive `REPRODUCED_DOWNLOAD_MANIFEST.json`; it no longer requires the missing external temporary metadata path or overwrites old manifests.

The T-slope accuracy diagnostic uses the source's stated0.35C instrument accuracy as a conditional envelope; source implementation comments reference0–50C. It is not a certified uniform bound for all environmental conditions. Neither this qualification nor the scalar-Ri diagnostic converts unresolved profile gradients into valid flux measurements.
