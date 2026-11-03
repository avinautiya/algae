# Inverse heat transport: development, identifiability and independent data gates

This package follows the failed all-baseline Pirámide field comparisons. It implements and checks stronger physical inference, rather than treating a temperature fit as proof of correct melt. **No new field superiority, Nature-level discovery or historical mathematical originality is claimed.** The second tower's camera outcomes remain unopened.

## Implemented work

- Conservative heterogeneous conduction and training-only objective comparisons in `physics/`.
- Exact warm-regime heat-scale ambiguity and sensitivity diagnostics, including finite-difference rank artefacts.
- Conditional cold-ice reference admittance inversion and interval parameter enclosures in `theory/ice_reference.py`.
- Positive integrated-energy conductance inference and a necessary conservative LP outer enclosure in `theory/weak_graph.py`.
- Primary-source novelty/identification review and independent transient-oracle tests in `literature/`.
- Verified Himalayan NgozumpaNG2 and Imja thermal observations and source-phase eligibility checks in `data/`.

## Actual development findings

Two-layer temperature fitting reduces raw TRAIN melt from371.8 to264.0mm water equivalent, still far above the already known85.2mm TRAIN observation. The fitted warm temperatures can remain identical while plausible capacity/conductivity scaling changes melt from170.6 to432.2mm. This is an exact model symmetry, not evidence that the field itself has that exact uncertainty interval.

The new integrated energy estimator recovers heterogeneous conductances on independently generated synthetic data. On real TRAIN profiles its inferred deepest conductivity changes from2.81 to0.419W/m/K when window length changes from6h to24h. With the stated conditional sensor bounds and assumed conductivity prior [.05,6]W/m/K, a necessary LP enclosure still permits approximately [.244,6]W/m/K on that edge. Its upper endpoint is a prior restriction, not a measured bound. The returned edge heat is not automatically heat delivered to the ice.

Cold-ice coupling can break the warm absolute-scale symmetry under strict reference assumptions. The provided inversion requires complete nonmelting history, known initial/deep ice conditions, independent ice properties and a known interface observation operator. It refuses to treat isolated cold snippets or unknown initial ice as sufficient. This is established admittance mathematics; the source audit identifies close prior work and no automatic historical novelty.

Pirámide's deepest TRAIN probe never falls below its stated±.4°C accuracy, so it supplies no resolved cold-interface reference. Actual NgozumpaNG2 data reach−6.18°C against stated±.25°C and provide approximately107days of fully cold lower debris. However, below-interface ice profiles and independent paired melt are absent; the conditions required for absolute cold-reference inference remain unresolved. Imja's strong cold spike occurs after the published field-period end and is not accepted as verified interface cooling.

## Reproduce

The sibling package `../piramide_thermal_test_20261103` is included in the repository branch. If using separate ZIPs, set `PIRAMIDE_RESEARCH_ROOT` to its extracted path. The dependencies are pinned there in `data/requirements_decode.lock`. Local private-temp paths are only an optional existing-workspace fallback.

Run `python3 verify.py`. `python3 theory/reproduce_training_graph.py NEW_OUTPUT.json` reproduces the integrated-energy training diagnostics and refuses to overwrite an existing output. Other study scripts also preserve their existing development outputs; use an isolated copy for a rerun. Portability-only path changes after calculation are recorded separately from numerical results.

## Research decision

This is development evidence and an executable inference foundation, not a successful new melt benchmark. Layering, goal-oriented inversion, thermal-admittance identification and RC-network regression have substantial prior art. The potential specific contribution is conservative inverse-flux feasibility with source-error controls and honest abstention, followed by a stable, independently blocked test. A new test should proceed only after the estimator's practical identification and claimed novelty survive the documented gates. The four previously scored periods must remain in the record.
