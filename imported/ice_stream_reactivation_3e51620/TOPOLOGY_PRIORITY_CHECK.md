# Coarse-cover topology and false contact tipping

Focused primary-source check, 31 October 2026. Claude work remains untouched.

## Candidate

At equal total coarse-debris cover fraction and solid volume, compare a longitudinal connected bare-ice corridor with a transverse spanning shield. A scalar fractional-cover closure assigns them the same thermal resistance, whereas a spatially resolved routing model may predict different incision and stream–cliff contact. The potentially useful contribution is a **connectivity-aware upscaling rule that prevents false contact-loss predictions**, with a verified regime of applicability.

No exact primary study of this equal-cover arrangement experiment, including stream reactivation and a demonstrated scalar-cover false tipping threshold, was found in targeted searches. That absence does not establish priority. A topological description alone is also not new mathematics: series/parallel resistance and flow connectivity are classical.

## Closest primary sources checked

### Mantelli, Camporeale & Ridolfi (2015)

[Supraglacial channel inception: Modeling and processes](https://doi.org/10.1002/2015WR017075). Primary full text §§2–3 and §8 inspected.

The model treats a distributed turbulent film over initially flat ice, with spatially uniform roughness, using a Stefan-like morphological stability analysis. It explains channel spacing through hydrodynamic feedback. Sediment is not a prognostic state; the conclusion explicitly leaves sediment–hydrodynamic deposition and localized melting as an extension. Thus it is a clear precursor for thermal channelization but does not supply a coarse-patch connectivity law or model the proposed arrangement counterexample. Do not present preferential thermal channel formation as newly discovered.

### Wytiahlowsky et al. (2025)

[Distribution and characteristics of supraglacial channels on mountain glaciers in Valais, Switzerland](https://doi.org/10.5194/tc-19-6461-2025). Primary HTML §§4.2–5.3 inspected.

The mapping study distinguishes bare, proximal patchy-debris and continuously covered channel settings. Channels occur along bare-ice/debris interfaces and moraine-adjacent depressions; patchy-debris channels are associated with greater sinuosity. These observations support the environmental relevance of arrangement and interfaces. The paper does not solve a topology-dependent routing or reactivation model. It also makes qualitative corridor self-selection unsurprising: the contribution would need to quantify predictive error and supply a useful closure rather than announce that water follows ice–debris gutters.

### Existing debris heterogeneity modeling

[Nicholson et al. (2018), Supraglacial debris thickness variability: impact on ablation and relation to terrain properties](https://doi.org/10.5194/tc-12-3719-2018); [Reid et al. (2012), Including debris cover effects in a distributed model of glacier ablation](https://doi.org/10.1029/2012JD017795). Primary abstracts/model framing inspected, not full equation-level reads in this focused check.

These are obstacles to any broad “first spatially heterogeneous debris” claim. They do not, on the evidence inspected, establish the exact open-water connectivity/reactivation experiment. Before claiming originality against either, their complete methods and relevant supplements must be checked.

### Hydraulic transmission beneath a shield

[Irvine-Fynn et al. (2017), Supraglacial ponds regulate runoff from Himalayan debris-covered glaciers](https://doi.org/10.1002/2017GL075398). Primary interpretation of through-debris pathways inspected.

The observations and interpretation include hydrological connectivity through thawing porous debris. Therefore a coarse transverse strip cannot automatically be represented as a perfectly blocking wall. Hydraulic transmission through pore space or an under-rock gap can retain water contact even where no open bare corridor exists. The problem is a coupled open-water/porous transmission problem unless the impermeable-shield limit is independently justified.

[Boyd, Goetz & Ham (2004), Supraglacial stream incision into debris-covered ice, Matanuska Glacier](https://studio.m-anage.com/gsa/2004NC/webprogram/Paper71537.html), primary conference abstract read. It documents sediment-bearing ice streams and lateral/vertical incision; it is not an equal-cover topology study. The original GSA endpoint failed; the conference-hosted mirror returned the abstract. Do not confuse a later bibliography's 2024 typo with the actual 2004 meeting date.

## Conditions needed for a bare-corridor comparison theorem

A lower bare corridor does not automatically remain dominant in every environment. A comparison argument needs all of the following declared or proven:

- matched route lengths, slopes, downstream boundary conditions and the same per-width hydraulic law;
- compatible width-dependent friction, especially sidewall effects in a narrow corridor;
- equal upstream water temperature, with cooling and atmospheric heat exchanges modeled consistently;
- monotonic integrated ice heat flux with per-width discharge under the chosen closure;
- geometry evolution that cannot introduce an upstream control or lateral cutoff absent from the reduced route model.

For a prescribed heat-transfer coefficient `h(q)` that is nondecreasing in per-width water flux `q`, a temperate-ice segment with no extra heat source has integrated heat removal proportional to

`q [1 − exp(−a h(q)/q)]`.

This expression is nondecreasing in `q`: increasing mass flux both supplies more sensible heat and reduces residence time, with the first effect sufficient under that monotone closure. This supports one restricted thermal comparison; it is not a proof for arbitrary Navier–Stokes/free-surface geometry.

## Decisive implementation test

Resolve transverse exchange and porous leakage before claiming that arrangement changes physical contact. Run the same actual spatial initial conditions through (i) fractional-cover averaged, (ii) connected-corridor reduced and (iii) independent spatial flow–heat reference models. Report errors in reactivation discharge, time to contact loss and cumulative heat delivered to the face. Include a regime where the averaged model works. A model that merely tags cells as blocked by design cannot independently establish the hypothesized connectivity threshold.
