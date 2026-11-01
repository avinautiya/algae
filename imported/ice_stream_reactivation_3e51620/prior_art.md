# Environmental mechanism priority audit

31 October 2026. Independent research; Claude files and running work are untouched.

## Decision

The defensible open question is **how much ice-cliff lifetime and melt change when a stream has enough heat to retreat the ice but insufficient capacity to remove externally supplied, coarse mantle debris**. The mechanism is environmentally motivated, but neither stream maintenance, shielding, transport limitation nor survival thresholds are new by themselves. A successful contribution needs a quantitatively different predictive system and a falsifiable result beyond transferring a permafrost model.

## Primary sources and exact obstacles

### Petersen, Hock & Loso (2024)

[Stream hydrology controls on ice cliff evolution and survival on debris-covered glaciers](https://esurf.copernicus.org/articles/12/727/2024/esurf-12-727-2024.html), DOI 10.5194/esurf-12-727-2024. Read primary HTML §§5.1–5.3 and methods/conclusions.

The paper already documents stream undercutting, debris clearing, avulsion-related cliff abandonment and subsequent burial. Section 5.3 explicitly proposes the exposed-cliff–runoff–stream-maintenance feedback, and blockage by larger clasts. It estimates smooth-ice clast mobilisation through force balance. Section 5.2 leaves the relative importance of sediment erosion/transport versus thermoerosional undercutting unquantified. Thus a first-hydrological-feedback claim is untenable; a predictive attribution and lifetime law remains a candidate gap. Observational association is not a calibrated dynamic transport model.

### Buri et al. (2016)

[A physically based 3-D model of ice cliff evolution over debris-covered glaciers](https://doi.org/10.1002/2016JF004039), DOI 10.1002/2016JF004039. Primary publisher methods/results accessible in search excerpts; not full equation-by-equation read here.

Already models evolving three-dimensional cliffs, atmospheric ablation, adjacent ponds and reburial. A geometric surrounding-debris view-factor threshold gates reburial. Therefore “first evolving cliff” or “first reburial model” would be false. A dynamic mass-conserving, size-resolved toe stock driven by stream removal should be compared against this geometric treatment rather than a straw-man fixed cliff.

### Evatt et al. (2025 publication; DOI dated 2024)

[A reassessment of ice cliff dynamics upon debris-covered glaciers](https://doi.org/10.1017/aog.2024.34). Read primary HTML model equations (1)–(17).

Already derives the moving-frame retreat relation, quasi-steady slope selection and critical debris thickness. It explicitly discusses unequal upper/lower sub-debris melt and changing cliff length. No stream forcing term was found in the inspected model. A new analysis must not rediscover moving-frame corrections or call ordinary slope selection new. It could instead determine whether a stream-maintained state is dynamically reachable, or is eliminated by debris delivery and hydraulic disconnection, which its quasi-steady construction does not resolve.

### Moore (2018)

[Stability of supraglacial debris](https://www.morageology.com/pubs/298.pdf), DOI 10.1002/esp.4244. Read primary PDF model-development portions, pp.287–288, plus introduction.

Already combines Darcy flow, ablation-generated meltwater and limit-equilibrium slope stability, deriving bounds depending on melt, permeability, thickness and slope. It assumes steady stresses/pressures and discusses transients as perturbations. Do not invent a moisture-dependent angle of repose and call that the first hydromechanical debris instability. A dynamic toe influx can use this as a supply closure, with assumptions stated.

### Douglas & Lamb (2024), plus Douglas et al. (2023)

[A model for thaw and erosion of permafrost riverbanks](https://doi.org/10.1029/2023JF007452); [author PDF](https://lamb.caltech.edu/documents/34126/JGR_Earth_Surface_-_2024_-_Douglas_-_A_Model_for_Thaw_and_Erosion_of_Permafrost_Riverbanks.pdf). Read primary publisher theory development and conclusions. Related [entrainment/slump-block paper](https://doi.org/10.1029/2023GL102974).

This is the strongest mathematical-mechanism obstacle: two moving fronts, a thawed insulating sediment layer, entrainment-limited/thaw-limited/buffered erosion, and thickness/failure cycles already exist. Their warming-sensitivity distinction precedes a glacier transplant. The material source differs: a largely clean ice cliff receives debris from an external supraglacial mantle, whereas thaw of sediment-bearing permafrost produces sediment locally. That distinction must cause a demonstrably new prediction; changing labels is insufficient. A scalar resistance-plus-clearing ODE alone is likely an application rather than a major discovery.

### Current glacier-scale debris models

[Verhaegen & Huybrechts (2026)](https://doi.org/10.1029/2025JF008748) already couples higher-order three-dimensional ice flow to englacial/supraglacial debris and mass-conserving margin evacuation. Read primary abstract/conclusions; not a complete source-code audit.

[Hardmeier et al. (2026)](https://tc.copernicus.org/articles/20/5365/2026/tc-20-5365-2026.html) evolves debris forcing and glacier response, while listing cliff formation/impact and more complete surface redistribution among future improvements. Read primary abstract/conclusions. Do not claim all glacier models omit debris transport; local stream-controlled cliff lifetime is the narrower gap.

## Candidate worth testing, inferred from the above

Use a **grain-size-selective, externally supplied toe stock coupled to stream hydraulic connectivity**. The important state is not simply blanket thickness: coarse grains may divert flow away from an ice face before the total debris volume is enough for geometric reburial. This permits a testable distinction between:

1. thermal protection at a still-connected toe;
2. hydrological disconnection by coarse load;
3. physical face burial.

Potential output: time to stream–cliff disconnection, time to loss of exposed ice, and cumulative excess melt. Compare the full system against transport-only armouring, thermal-only erosion, geometric reburial and a matched Douglas–Lamb-style scalar buffer. Do not assume the full system wins.

A novel environmental consequence would require showing that size distribution or forcing order changes one of those outputs substantially at fixed total debris supply, water and heat budgets, under independently resolved hydraulics. The specific consequence is **a hypothesis to compute**, not an established discovery or verified global gap. Coarse-grain diversion already appears qualitatively in Petersen; novelty must be quantitative prediction and mechanism discrimination.

## Necessary no-shortcut tests

- Keep debris as solid volume; distinguish mantle influx, newly melted-out grains, export and storage. Avoid double-counting mantle debris as an ice-volume concentration.
- Use an actual local flow/grain obstruction solution or independently calibrated closure to predict connectivity. An arbitrary logistic switch can manufacture bistability.
- Report hydraulic work/transport capacity and thermal energy separately. Heat flux and discharge are not interchangeable.
- Compare short pulses and long steady forcing with identical water/heat/sediment totals. A response difference must survive temporal/spatial refinement.
- Use the same observation operator for reference and baseline: exposed true area, outline area, stake-normal ablation and horizontal retreat differ.
- A simulated closure correction is computational validation within scope; it is not proof of field skill without an independent environmental dataset.

## Follow-up: Exner–Stefan routing rather than a prescribed abandonment switch

A stronger candidate supplied by the parent is a toe channel and an ice-floored bypass, coupled through upstream water storage. Each branch has ice elevation `z_i`, sediment solid stock per plan area `s_i`, porosity `p_i`, and hydraulic bed `eta_i=z_i+s_i/(1-p_i)`. Resolve water partition from hydraulic heads, not a logistic switching function. Then

`d eta_i/dt = -E_ice,i + (D_i - div q_s,i)/(1-p_i)`.

Here `D_i` is solid-volume deposition per plan area, including external mantle supply associated with cliff retreat. Thermal incision is obtained from the actual heat flux to ice after sediment resistance and water cooling. Raising the toe bed reduces its water share, which can reduce both clearing and ice incision; the lower bypass may then capture the stream. This equation alone is a conservation identity, not novel mathematics. Predictive content is in a verified thermal/hydraulic/sediment closure and its environmental consequences.

### Additional primary obstacles

[Bolla Pittaluga, Repetto & Tubino (2003), Channel bifurcation in braided rivers](https://doi.org/10.1029/2001WR001112), inspected primary publisher formulation §§2–3 and stability §6. Already couples water discharge, sediment continuity and transport capacity at bifurcations, finding asymmetric equilibria and branch closure. A water-split/Exner instability is not a new mathematical invention.

[Kleinhans et al. (2008), Bifurcation dynamics and avulsion duration in meandering rivers by one-dimensional and three-dimensional models](https://doi.org/10.1029/2007WR005912), primary abstract and methods excerpts inspected. Explicitly builds on bifurcation sediment partitioning and evolving branch dominance. Full paper not read here; do not make an exhaustive exclusion claim.

[Karlstrom et al. (2013), Meander formation in supraglacial streams](https://doi.org/10.1002/jgrf.20135), primary abstract/model framing inspected. Already models thermal surface/channel erosion in supraglacial stream migration, while assuming little sediment and neglecting relative channel incision in that study. The ice-substrate component is therefore not new either.

### Potentially differentiating test

At fixed upstream water and heat inputs, compare the coupled ice-bed/debris model with (a) frozen ice substrate + Exner, (b) thermal ice + no sediment stock, and (c) externally prescribed stream disconnection. A consequential candidate result is that **mantle supply generated by cliff retreat destabilizes a stream-maintained cliff even while direct thermal melting increases**, because route disconnection removes its maintenance flow. The critical test is whether this occurs with physically independently justified supply and transport rates, persists with resolved hydraulic routing, and differs from known alluvial avulsion and permafrost shielding limits. No exact prior supraglacial model of this full combination was found in the targeted search; this is not proof of first occurrence.

The environmental outputs should be quantitative: critical cumulative delivered coarse volume, disconnection/rescue discharge threshold, lead time to route loss, and lost cumulative excess cliff melt. A figure showing only phase-space features without these predictions would not meet the user's objective.
