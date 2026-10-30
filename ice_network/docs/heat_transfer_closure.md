# Heat-transfer closure at the melting wall: what is assumed, what is known, what the results depend on

**Status: verified computational model with an UNVALIDATED open-channel heat-transfer closure.** No experiment was obtained or required for this milestone; no claim below is a physical validation.

## What the code uses
h_T = Nu k_w / D_h, Nu = nu_factor x 0.023 Re^0.8 Pr^0.3 (Dittus-Boelter with the exponent for a fluid being cooled), D_h = 4A/P, Re = v D_h / nu, water at 0 C (Pr 13.4). Heat flux is uniform over the wetted perimeter. The same expression is used (separately coded) by the resolved reference solver, on purpose, so that the reduced-versus-reference comparison tests the reduction and not the closure.

## Regime in which it is a defensible starting point
Fully developed, turbulent, constant-property flow in a smooth conduit, nominal 1e4 <~ Re <~ 1.2e5, 0.7 < Pr < 160. In the controlled cases here Re ranges from about 2e4 to 1.7e5 (`diagnostics.applicability`), so the upper end is outside the nominal pipe range. **Not represented:** free surface and secondary currents of an open channel, entrance (developing-flow) effects over short reaches, ice roughness and scallops (which change both friction and heat transfer), a wall whose geometry and temperature are being changed by the very heat flux, variable properties near 0 C.

## Alternatives that exist, and how far they differ
Established smooth-pipe correlations already disagree at these Reynolds numbers. Computed here at Pr = 13.4 (`records` not needed, closed form):

| Re | Dittus-Boelter n = 0.3 (used) | Dittus-Boelter n = 0.4 | Gnielinski (Petukhov friction) | Gnielinski / used |
|---|---|---|---|---|
| 1e4 | 79 | 103 | 101 | 1.27 |
| 2.3e4 | 155 | 200 | 216 | 1.40 |
| 5e4 | 288 | 373 | 428 | 1.49 |
| 1e5 | 501 | 649 | 786 | 1.57 |
| 1.7e5 | 766 | 993 | 1252 | 1.63 |

So the choice of pipe correlation alone spans a factor of about 1.3 to 1.8 in the transfer coefficient. Roughness is expected to push transfer higher still, and a free surface can reduce it; the factor 4 range of `nu_factor` (0.25-4) used in the sensitivity study is a bracketing device, not a literature value. Open-channel or ice-wall specific correlations (e.g. from experiments on melting ice surfaces in flowing water, or friction-velocity based Stanton-number closures used in glacier-channel work) were **not located or read in this milestone**; identifying them and the geometry and Reynolds range each was measured over is the next literature task and may allow narrowing the band.

## How much the predictions depend on it (computed, `records/incision_experiment.json`, E3)
Warm inflow (1 C, 0.05 m3/s), 100 m reach, 8 h, converged numerics (stations 17, dt 450 s):

| nu_factor | sill lowering (m) | heat to ice (GJ) | delivered volume (m3) |
|---|---|---|---|
| 0.25 | 0.0144 | 0.363 | 979.7 |
| 0.5 | 0.0288 | 0.701 | 1004.2 |
| 1 | 0.0572 | 1.311 | 1051.1 |
| 2 | 0.1116 | 2.322 | 1135.2 |
| 4 | 0.2070 | 3.748 | 1261.2 |

Incision is almost proportional to the multiplier (exponent about 0.95 at low values, saturating at high exposure). A factor-2 closure uncertainty therefore means a factor-2 uncertainty in predicted incision and a 4-8 % change in delivered volume in this case, while the numerical errors of the converged scheme are about 1 %. **The closure, not the numerics, dominates the uncertainty of the warm-water incision prediction.** This is a sensitivity result, not validation.

## What would be needed to narrow it (without field campaigns)
1. Published laboratory measurements of heat transfer to a melting ice wall under turbulent open-channel or duct flow, with geometry, Re, temperature and roughness stated; compare against the used correlation over their range. Existing experiments may serve as optional validation; they are not required for the core project.
2. A resolved calculation (RANS or LES with a heat-transfer wall treatment) of a straight, fixed-geometry, warm-water channel with a melting-point wall, to calibrate an effective multiplier and its Re-dependence; this was scoped, not launched (no CFD campaign in this milestone).
3. Until then the results are labelled conditional on the closure, with the band 0.5-2 used for reported sensitivities.

## Manning roughness
The Manning n = 0.04 is assumed, shared by both models, and is an input uncertainty of the same kind: it enters the discharge for a given stage and the friction dissipation. Its effect is quantified in `docs/evaluation_results.md` where run, and is separate from the heat-transfer closure.
