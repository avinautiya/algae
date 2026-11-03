# Signed basal control-volume development candidate

**Verdict: it fails the useful-improvement gate.** The physically motivated storage correction is too small to explain the released prediction errors. It also produces negative half-hour interface-energy increments, so the signed net energy must not be presented as a verified no-refreezing melt trajectory. No Tower3 or later Tower1 outcome was read. T1/T2 results are now DEVELOPMENT only.

## Equations

Take a static basal control volume from the midpoint between DT4 and DT5 down to a temperate0C ice boundary. Let h be their spacing and g the positive assumed DT5-to-ice gap. Piecewise-linear temperatures give stored sensible energy per area

    H = C [h*T4/8 + (3h/8 + g/2)*T5].

The midpoint temperature is (T4+T5)/2. Independent trapezoidal integration of its two linear subsegments gives the same expression. The top conductive flux approximation is k(T4−T5)/h. With no advective, latent-pore, lateral or volumetric source,

    Q_ice(I) = k integral_I (T4−T5)/h dt − [H(end)−H(start)].

This is an integrated control-volume energy identity under explicit assumptions, not new mathematics. Signed temperatures and fluxes are retained. Q/L has units kg/m², numerically mm w.e., with L333700J/kg. The per-site conductivity is calibrated only from that site's already-known TRAIN melt:

    k = [L*Y_train + DeltaH_train]/J_train.

No development-evaluation target is used to fit k. All period endpoints match the released native camera operators (including Tower1's previously amended Dec30T00 endpoint). The five native14-day intervals use complete half-hour temperature support; gradient integration uses signed trapezoidal interpolation.

## Development results

| Candidate/baseline | T1 four-period MAE | T2 four-period MAE |
|---|---:|---:|
| Signed gradient alone |16.712128|7.716660|
| Storage correction C1.855487MJ/m³K,gap0.02m |16.478587|7.675989|
| Air PDD (released comparison) |7.906105|7.265305|
| Calendar-modulated PDD (T1 development fit; T2 transfer) |2.252768|3.434992|

The primary candidate's inferred conductivities are1.511605W/m/K at T1 and3.020456W/m/K at T2. This large site dependence is a conditional calibration result; it is not a transferred material law or proof of true conductivity. Local gradient operators, unresolved position, nonconduction and camera calibration all enter the effective value.

All24 development trials are retained: capacities1.5,1.855487,2.5MJ/m³K and basal gaps0,.01,.02,.04m, at both sites. The capacities and gaps are assumed finite sensitivity cases, not measured uncertainty distributions or a certified source envelope; the T1-derived central C proxy is not a measured T2 C. The zero-gap row is a limiting calculation. Source depths and thickness carry separate uncertainties. No winning case is promoted: even the lowest trial MAEs,16.368508 at T1 and7.653354 at T2, fail air/calendar comparison.

At the primary values, endpoint sensible-storage changes are only−0.440657 to+0.319525mm w.e. at T1 and−0.188760 to+0.139071mm w.e. at T2. Corrections this small cannot fix10–30mm errors. Increasing C to manufacture a large correction would contradict its independently bounded physical role.

## Causal phase-conversion gate

The fine-step ledger gives259 negative interface-energy increments of3360 at T1,226 of3360 at T2. Their aggregate signed energies correspond to−3.389781 and−3.301214mm w.e. These were not clipped away. Positive full-period net Q does not prove all supplied energy went into melt at each time step.

A genuine causal melt model must allocate negative heat to an explicitly modeled ice cold-content state, or justify available water/refreezing and its mass ledger. Simply taking max(Q,0) creates additional melt without the compensating cold/storage history. A reflected cumulative-energy reservoir can track cold content algebraically, but its basal temperature and coupling back into the unmeasured gap must then be derived consistently; this candidate's fixed0C P1 boundary does not already supply that derivation. Therefore the reported development numbers are signed conductive melt-equivalent diagnostics, not a fully physical phase-resolved field predictor.

## Tests and files

Six meaningful tests pass in `test_storage_flux.py`: independent P1 storage integration, negative signed-gradient preservation, exact energy ledger/TRAIN recovery, release of stored heat during cooling, explicit negative-energy flags, and invalid-input refusal.

`storage_flux.py`: executable core.
`run_development.py`: complete bounded development calculation; deliberately uses exclusive output creation.
`DEVELOPMENT_STORAGE_RESULTS.json`: all24 trials, no winner filtering.
`DEVELOPMENT_STORAGE_PHASE_DIAGNOSTIC.json`: original results plus fine-step phase gate.

Reproduce the core tests with `python3 -B -m unittest discover -s algorithm -p test_storage_flux.py`. Existing outputs are never overwritten by the runner. Runtime is under1second; no heavy computation or new data download was used.

## Contribution and next decision

The useful result is a quantitative falsifier: endpoint sensible storage in this basal control volume is too small to rescue the current field errors under defensible capacities. This prevents an unsupported expensive campaign or artificially large fitted storage coefficient. It is not a new environmental discovery, validation success or first use of heat storage in glacier debris. Temperature-only phase/transport closure remains unresolved; the failed field models and all frozen results are preserved.
