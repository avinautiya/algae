# Rochemelon 2005 external observation-match gate

## Decision

**Not ready for an unqualified independent field-accuracy score.** A conditional breach-incision comparison may be constructed, but section identity, survey timing and hydraulic assumptions must first be explicit. No numerical bed-target array or new model score has been extracted in this audit. Root's external parameter freeze preceded visual inspection of the target-bearing figure.

## Primary observations and source identity

Vincent et al. (2010), DOI [10.3189/002214310791190857](https://doi.org/10.3189/002214310791190857), describe repeated theodolite/staff bed surveys at sections A-C and less frequent downstream profiles. Figure 6b's black dots are measurements, whereas its continuous curve is lake level. Figure 5 identifies spatial sections and dated profiles; exact survey hours are not supplied there. Lake thermistors have ±0.1°C accuracy. Discharge is derived from lake-level change, hypsography and recharge, with reported ±0.1 m³/s accuracy. Generic geodetic precision is not a separately stated error for every bed target. Cold ice, mechanical ablation and dissipative energy enter the study context. The authors fitted their heat-transfer coefficients to this event, so those coefficients cannot supply an independent event benchmark.

Source PDF SHA256 and visual inspection pages are recorded in `OBSERVATION_MATCH_GATE.json`. Methods page 93 and figure-design page 94 were inspected at readable resolution. The existing full source remains in the sibling forcing folder; it is not edited.

## Binding and scoring requirements derived for this project

1. **Bind the target to a location.** Figure 6b has one black-dot bed series without A/B/C labels in its caption or visible legend. The paper's breach model suggests an upstream interpretation, but that is not a demonstrated station binding. Do not duplicate it into three apparent independent station series. A conditional breach-point interpretation must be labelled as an assumption unless a primary source resolves it.
2. **Use an observation-equivalent output.** Positive incision is initial vertical bed elevation minus later elevation at the same fixed section. A plan-view volume, a length-averaged lowering, wall recession and total ice melt are different observables. They cannot be compared directly to these point elevations. Mechanical removal, if present, is not separately measured thermal melt.
3. **Keep time support honest.** The common elapsed-hour figure origin avoids inventing a timezone. Calendar-dated spatial profiles do not establish exact-hour forcing windows. If their times remain unknown, evaluate a declared timing envelope; do not assign precise noon timestamps and count that assignment as observed timing.
4. **Fix initial conditioning before target extraction.** If the first eligible observed bed point initializes the simulation, exclude that conditioning point from accuracy scoring and retain every later eligible point. Do not fit a vertical offset to all target points. A construction-based post-initial-period start rule may be justified before error inspection, but the same rule must apply to all candidate models.
5. **Respect the inlet/local distinction.** Lake temperature is an inlet thermal boundary. Applying it unchanged at downstream B or C requires a cooling/transport assumption. Likewise, width/depth inferred from critical-flow hydraulics remain inferred inputs. A model supplied with measured future stage/geometry tests a conditional melt closure, not its ability to forecast those variables.
6. **Separate energy and geometry assumptions.** Transferring a local Swiss convective closure alone does not account for gravitational dissipation, cold-ice heat uptake, changing cross-sections or spatial water cooling. Their neglect or approximation must be declared before scoring. The source's event-fitted coefficients may be shown only as an in-sample historical reference, never as independent predictive skill.
7. **Keep uncertainty conditional.** Pixel/axis uncertainty can be propagated after digitization. It does not replace an unreported target-measurement error. Instrument bounds and deterministic scenario envelopes are not statistical confidence intervals. Repeated epochs at one drainage event are not independent glaciers.

## Release gate

Before extracting targets, record the selected section interpretation, allowable timing envelope, initial-state rule, common endpoints/footprint, parameter freeze, treatment of inlet cooling/dissipation/cold ice, and source-error assumptions. If section binding cannot be resolved, a clearly labelled conditional sensitivity study is possible; a precisely matched external validation claim remains unsupported. Missing raw tabular observations must not be invented from smooth model curves.
