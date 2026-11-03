# Finite radiation/cooling reservoir development study

All Tower1andTower2targets used here were already released. This is DEVELOPMENT ONLY, including every timescale trial and retrospective resampling diagnostic. No Tower3or later Tower1camera data were accessed. Previous frozen results remain unchanged.

## Effective model and inputs

For fixed source-calibration values albedo a=.2and emissivity epsilon=.94, define

Q_proxy=(1−a)SWIN+epsilon*[LWIN−sigma*(T_buried+273.15)^4].

Incoming radiation comes from the previously available ERA5-labelled weather input. The native near-surfaceDT1thermistor is used as an outgoing-longwave temperature proxy. It is buried by approximately2cm; it is **not a measured radiometric surface temperature**. Turbulent sensible/latent heat, surface wetness, snow and the actual ice thermal profile are missing. Q_proxyis therefore incomplete atmospheric forcing, not observed net heat into the glacier.

An effective reservoir obeys

E'=Q_proxy−D_tau(E), D_tau(E)=max(E,0)/tau.

PositiveErepresents effective stored energy and negativeEan effective deficit. Negative forcing is retained; later positive forcing repays a deficit before energy is delivered. Neither sign is asserted to be measured glacier-ice cold content. For positiveEthe response is a stable exponential filter. The negative branch is an integrator, so the full system is not claimed globallyBIBOstable under indefinitely negative forcing. Only its finite source-forced response is studied.

InitialEuses source-profile temperature integration and the existing training-only heat-capacity proxy. Tower2's unobserved final1cm is extended with the last debris temperature for this initialization, not imposed as0°Cice. The run starts24November, before the2Decembertraining interval. Timescales were declared before this campaign as0,6hours,1day and3days at0.5m thickness; each site scales tauwith(d/.5)^2. All four trials are retained, including the zero-timescale limit.

## Exact piecewise-constant update and verification

WhileE>0under constantQ,

E(t+s)=Q*tau+[E(t)−Q*tau]*exp(−s/tau),
G(t,t+s)=Q*s−[E(t+s)−E(t)].

IfQ<0the positive branch hits zero aftertau*log[1+E/(−Q*tau)]. Remaining negative forcing decreasesElinearly with no delivered heat. IfE<0andQ>0the deficit is repaid after−E/Q, followed by the positive solution. ForQ<=0andE<=0there is no delivery. The tau=0limit releases positive initial storage and uses an exact deficit-repayment rule; it is not simple clipping of negative forcing.

Each step satisfies DeltaE+G=Q*Delta t. Seven tests cover steady state, cooling/reheating repayment, zero-crossing events, exact semigroup composition, radiative equilibrium, the zero-timescale limit and invalid inputs. Maximum absolute per-step ledger residual in actual runs is2.33e−10J/m². This verifies the effective implementation, not the completeness of environmental heat accounting.

Source forcing is evaluated at30minute midpoints with native temperature/radiation interpolation; net radiation is constant within each update. Repeating at15minutes changes issued interval predictions by at most0.00868mm water equivalent, far below the development errors. Thus this campaign's poor performance is not explained by its integration step.

## Calibrated output and physical interpretation

Potential delivery inmmwater equivalent is sum(G)/Lf, withLf=334000J/kg. One training-only site gainA=Y_train/H_trainis required to lie in[0,1]; every trial passes, withAapproximately.110–.118. Future model output isA*H.

The remaining(1−A)*Gis an unresolved partition. The identity DeltaE+A*G+(1−A)*G=Q*Delta tis closed **by construction**. It does not independently establish that the unresolved87–89%is atmospheric loss or identify any missing physical mechanism. It must not be presented as observed energy closure or a measured latent fraction. A complete surface-energy model would need actual surface temperature, surface wetness and independently justified turbulent/phase terms.

## Results, all trials preserved

| Reference tau | Tower1periodMAE | Tower2periodMAE |
|---|---:|---:|
| 0days | 17.5952 | 13.6301 |
| .25days | 17.6048 | 13.6539 |
| 1day | 17.2562 | 13.4561 |
| 3days | 17.4143 | 13.3932 |
| Locked airPDD | 7.9061 | 7.2653 |
| Locked calendar | 2.2528 | 3.4350 |

Units aremmwater equivalent. The locked calendar parameter was originally selected using all fourTower1development observations, so that particular Tower1row is a fitted development score, not leave-one-period-out evidence.

For additional retrospective diagnostics, each evaluation period is left out of timescale selection over the other seven already-inspected periods; one local training gain remains allowed per site. The calendar baseline is also refitted under the same period exclusion using its previously declared parameter grid. PooledperiodMAEs are15.7205for the reservoir,3.6632for calendar and7.5857for locked air. Periods share observation endpoints and sites share glacier/weather: these are not independent replicate validation samples or statistical confidence intervals.

With timescale selected on the other stake's four periods, reservoirMAEis17.4143atTower1and13.4561atTower2, againstair7.9061/7.2653andother-stake-calendar2.3331/3.4350. These are retrospective transfer diagnostics, not new untouched field tests.

The reservoir changes the late-period bias direction but does not improve prediction. Its output varies little across the four timescales, so this finite family does not support the claim that adding one effective radiation/storage timescale fixes the errors. This is not a universal rejection of heat storage or radiation physics.

## What remains unestablished

The model uses an incomplete and partly reconstructed forcing proxy. A buried thermometer cannot independently supply true outgoinglongwave; deeper moisture probes cannot define surface evaporative availability. Missing snow values do not prove snow absence. Known full debris models already include atmospheric exchange and storage, so the present one-reservoir reduction has no historical novelty claim. No result here establishes a superior validated melt model or resolves the earlier failed comparisons.

The next justified work is a **new explicitly specified** atmospheric/ice boundary model supported by independent inputs, or a different observable such as the native cold-response temperatures. Do not silently add fitted terms to this failed family and reuse these same development errors as validation.

## Files and commands

`radiation_reservoir.py`, `test_radiation_reservoir.py`, `run_development.py`; all trial predictions, gains, refinement checks, source hashes and retrospective diagnostics are in `results_v1/DEVELOPMENT_RESULTS.json`. Eight source-forced traces preserve states and ledgers. Figure scripts render the complete family comparison and an explicitly illustrative effective-state trajectory.

```
PYTHONDONTWRITEBYTECODE=1 /private/tmp/glacier-round2-20261103/data/venv312/bin/python -m unittest discover -s alternative -p 'test_*.py' -v
PYTHONDONTWRITEBYTECODE=1 /private/tmp/glacier-round2-20261103/data/venv312/bin/python alternative/run_development.py alternative/results_new
```

The current executable uses the sibling development packages at their recorded temporary paths; the result records their hashes. It never opens any native camera MATfile and refuses overwriting an output directory.

## Subsequent forcing-source audit note

The parallel source audit found per-siteMATforcing with explicitDeltaGMT=−4and a midnight first date, while the existingCSVbegins at01:00and has empty timezone metadata. Exact row/time-bin correspondence is being checked independently. This study used sourceCSVlabels unchanged and assumed a common clock with process temperatures; that is conditional. No alignment was selected by melt-error fit, and these development results do not resolve the source clock. A future fullSEBtest must bind that clock before scoring, together with true surface-temperature and evaporation assumptions.

### Resolved source-label correspondence

The subsequent source-agent audit resolved the apparent start-label discrepancy: the MATLAB weather export has one extra initial midnight row. Across 4,187 identical naive labels, MATLAB pressure minus the provided CSV pressure is a constant site-elevation adjustment (standard deviation approximately 1e-8 hPa); shifting by one hour destroys that identity. The author's driver also explicitly matches equal labels. This supports equal-naive-label weather binding, while the physical timezone remains conditional. It does **not** demonstrate a one-hour phase error, and no alignment was selected from melt errors. The separate neutral-profile energy hypothesis and its source-only no-go gate are in `neutral_flux/PROTOCOL_AND_NO_GO.md`.
