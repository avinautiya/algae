# Adversarial novelty audit of stream–cliff recovery

Audit date: 1 November 2026. This is an independent literature and experiment-design contribution; no active Claude file, frozen result or process was changed. Source inspection ranges are stated below. Search absence is not proof of priority.

## Decision

The strongest defensible target is **a quantitative boundary between debris clearance and sustained thermal reconnection after stream abandonment**, rather than discovery that streams maintain cliffs, debris buries them, or thermal channels retain history. Those larger mechanisms already exist in the literature.

A precise research question is:

> After a sediment obstruction is removed, when can ordinary meltwater flow recover contact with an ice-cliff toe, and when has differential incision of an alternative route made removal insufficient?

The useful output is a predicted recovery threshold, finite recovery deadline, or required water/heat budget, checked against an independent spatial solver. The contribution would be strongest if a model using debris amount alone predicts recovery while a geometry-and-heat-aware model correctly predicts failure across withheld layouts and forcing histories. No current result establishes that environmental improvement, and the existing failed partition gate must be resolved before testing it.

## Closest primary sources and exact claim boundaries

### Kneib et al. 2022 — direct Himalayan precedent

[Sub-seasonal variability of supraglacial ice cliff melt rates and associated processes from time-lapse photogrammetry](https://tc.copernicus.org/articles/16/4701/2022/), especially Discussion §§5.2–5.3, inspected in full HTML for this audit.

The authors explicitly describe stream disconnection at the 24K cliff as switching off sediment evacuation, causing partial reburial. They also connect pond drainage with loss of basal thermoerosion. They recommend integration of debris redistribution, hydrology and cliff energy balance. Thus **stream disconnection → lost debris evacuation → reburial is not new**. The inspected discussion does not provide a post-clearance reattachment threshold, rescue-energy certificate or comparison of matched debris removal at different incision states.

### Petersen, Hock and Loso 2024 — direct stream/cliff precedent

[Stream hydrology controls on ice cliff evolution and survival on debris-covered glaciers](https://esurf.copernicus.org/articles/12/727/2024/), §§4–6 inspected in primary HTML.

The study reports stream-associated cliff maintenance, clast transport, avulsion, subsequent reburial and a conceptual geomorphic model. It recognizes that large clasts can obstruct channels and that exceptional flows can reorganize drainage. Crucially, the relative contribution of debris evacuation versus thermoerosional undercutting is left unresolved. **Neither stream-maintained survival nor flood mobilization can be the headline discovery.** A controlled quantitative separation of these two mechanisms is a defensible target, although a plan-view bed model alone cannot predict full cliff lifetime.

### Moore — permeability and coarse-toe extinction already proposed

[Stability of supraglacial debris](https://www.morageology.com/pubs/298.pdf), abstract and concluding ice-cliff implications inspected. The PDF records acceptance in 2017; use the final journal record for publication-year metadata before manuscript submission.

This develops Darcy/meltwater and slope-stability theory. Its cliff discussion explicitly contrasts coarse permeable debris that can accumulate at the foot and extinguish backwasting with finer debris whose retained pore pressure can sustain evacuation. **Coarse toe accumulation and permeability-controlled persistence are established hypotheses.** The new study must not rediscover that contrast under an Ergun-drag label. The open-channel routing and evolving alternative-route recovery problem is more specific than the inspected formulation.

### Buri et al. 2016 — dynamic cliff modeling already exists

[A physically based 3-D model of ice cliff evolution over debris-covered glaciers](https://doi.org/10.1002/2016JF004039), publisher abstract/model summary inspected; not a new equation-level audit here.

The model already evolves cliff geometry, reburial and pond influence and is compared against observed evolution. A broad first dynamic cliff model or first incorporation of water/reburial claim is invalid. The present proposed advantage is an explicitly resolved stream-routing recovery mechanism, not replacement of atmospheric/radiative cliff physics. Full methods and supplement remain necessary before asserting an absent equation or closure.

### Evatt et al. 2024 — geometry/frame theory already substantial

[A reassessment of ice cliff dynamics upon debris-covered glaciers](https://doi.org/10.1017/aog.2024.34), open primary PDF abstract and model framing inspected.

The paper uses a moving reference frame relative to lowering debris-covered ice and predicts slope selection and critical debris thickness using Baltoro parameters. Therefore moving-frame corrections, endogenous slope selection and a generic critical-thickness threshold are not new. Our fixed plan-view surface must not turn local hydraulic contact into cliff retreat or ice-volume loss without the relevant geometric observation operator and frame conversion.

### Delaney et al. 2025 — thermal hydraulic memory and sediment hysteresis

[Sediment transport capacity response to variations in water discharge in pressurized subglacial channels](https://tc.copernicus.org/articles/19/2779/2025/), §§1, 4–5 inspected in primary HTML.

Evolving melt-opened/creep-closed channel size produces hysteresis between discharge and sediment-transport capacity. Their primary setting is pressurized subglacial channels. This defeats claims that evolving ice geometry creates sediment-transport memory for the first time. A surface toe/bypass system has a different control geometry and no overburden creep closure, so a tested **post-clearance reconnection limit** could still be a distinct result. It must survive comparison with a model that already carries geometric state.

### Mantelli et al. 2015 and general pulse transport

[Supraglacial channel inception: Modeling and processes](https://doi.org/10.1002/2015WR017075) was previously inspected at equation level in `research/environmental_mechanism_2026_10_31/TOPOLOGY_PRIORITY_CHECK.md`; that prior audit is relied upon, not newly reread here. Thermal channelization and instability are precursors, not new results in this project.

[Humphries et al. 2012](https://doi.org/10.1029/2011WR010419) was previously inspected in the companion `PULSE_PRIORITY_CHECK.md`. Equal total water supplied in different hydrographs already produces different sediment responses. A same-totals/different-outcomes plot alone does not establish distinctive glacier physics.

## What could still produce a consequential result

The strongest candidate is a **clearance–reconnection mismatch with a moving thermal barrier**:

1. A deposit diverts water from a toe route into an alternative ice route.
2. The alternative route incises during abandonment; the toe may incise more slowly under debris.
3. A later flood removes debris, but water returns to the now-lower alternative route when ordinary flow resumes.
4. A static-bed counterfactual or a model retaining only remaining debris incorrectly predicts recovery.

This is a hypothesis; the chain must be demonstrated without imposing blockage or permanent route walls by construction. A verified calculation showing merely temporary flood rewetting is not recovery.

The environmental significance is that sediment export or optical exposure may fail as a proxy for re-established basal stream heat delivery. The practical forecast is a conditional recovery class and heat-delivery error, not an unvalidated glacier-wide melt improvement. Applying this to Himalayan climates needs realistic forcing ranges and geometries; Alaska observations establish process relevance, not Himalayan calibration.

## Minimal decisive comparison after Claude's verification gate

Keep initial ice geometry, actual grain solid volume, inlet histories and total heat fixed across models. Carry hydraulic storage, latent exchange, remapping, export and actual event time in the ledger. Compare:

- evolving ice and prescribed instantaneous removal;
- the same removal with differential incision disabled;
- the same post-removal geometry but debris-free from restart;
- porous transmission and open under-rock/around-rock routes;
- a geometry-state reduced predictor versus a debris-only predictor.

This separates three explanations: remaining material blocks recovery; irreversible substrate geometry blocks recovery; or numerical wetting/drying and boundary artifacts falsely block it. Include a case that recovers, a case that never loses contact, and a no-debris baseline. Do not select only dramatic abandonment cases.

Use ordinary-flow contact persisting for a predeclared duration as primary outcome, assessed after the flood ends. Report time-to-loss, time-to-sustained-recovery, cumulative **local** ice heat delivery and actual debris export. Threshold definitions are operational choices, so show sensitivity and event-time brackets. Grid and morphology-step uncertainty must be smaller than any model-performance claim.

## Red flags that invalidate the headline

- The effect disappears when the entrance/reference mismatch is corrected.
- The effect requires an impermeable strip even though realistic porous or lateral routes keep the toe wet.
- Fixed inlet/outlet sills become controlling artifacts as ice lowers.
- Removing debris restores ordinary contact once event scheduling is exact.
- A geometry-state baseline performs as well as the proposed predictor with less complexity.
- The result is a classical avulsion response of the frozen-bed model, with negligible additional thermal contribution.
- Claimed cliff survival changes are inferred from a plan-view floor model that never evolves the cliff face or toe-debris supply.

## Claim ladder

**Available now:** a literature-grounded hypothesis and falsification plan, plus prior conditional reduced-model demonstrations.

**After a successful verified spatial campaign:** a computational mechanism result and a tested conditional forecasting rule over declared regimes.

**Not supplied by that alone:** historical priority proof, validated real-glacier improvement, universal control of cliff survival, or a Nature/publication/competition guarantee.

The literature does leave a credible narrow modeling question. Its significance must come from a quantitatively important and robust answer, rather than from calling established constituent processes novel.
