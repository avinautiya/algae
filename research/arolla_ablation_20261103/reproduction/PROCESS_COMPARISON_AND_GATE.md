# Process comparison and a decisive independent test

This is a source audit plus a target-free calculation. It establishes neither a new melt law nor field superiority. No held-out melt prediction was computed here, and no current model was retuned.

## What the professional model actually represents

The [final DCG-MIP supplement §2.1, PDF pp3–5](https://tc.copernicus.org/articles/20/1895/2026/tc-20-1895-2026-supplement.pdf) specifies separate snow and ice processes for DEBCF, a ten-layer debris conduction model and an iterative surface balance. Its evaporation activation uses precipitation, and its turbulent treatment includes Richardson-number limits. The archived matrix's exact exported snow/ice component remains unresolved. Separate equations in a paper do not identify the variable written to a particular CSV. Do not assert contamination or certify a like-for-like professional skill score until that mapping is authoritative.

The original [Reid and Brock 2010 open author PDF](https://www.pure.ed.ac.uk/ws/portalfiles/portal/11434574/Energy.pdf) offers MATLAB code by request. I did not locate the exact DCG-MIP executable version. The Zenodo outputs archive central directory contains matrices and date files, not model/analysis source. No author was contacted.

**Novelty exclusion:** adding longwave radiation, latent heat, a moisture state, or stability corrections is not itself new. [DCG-MIP §6.2](https://tc.copernicus.org/articles/20/1895/2026/) already compares these treatments and calls for direct turbulent-flux validation. [Steiner et al. 2018](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2018.00144/full) measures evaporation over Himalayan debris. A new contribution needs a demonstrably different prediction, not another implementation of this established balance.

## Actual target-free environmental contrast

`process_diagnostics.py` reads only Pirámide forcing from the original local ZIP, hashes archive and CSV, refuses missing support/sentinels, and writes `PIRAMIDE_PROCESS_DIAGNOSTICS.json`. The windows are the source-label dates 2014-12-04–2015-02-10 and 2015-02-10–2015-04-08. These are not claims about unknown stake visit hours.

| Quantity | First window | Second window |
|---|---:|---:|
|Complete hourly rows|1632|1368|
|Mean net shortwave W/m²|276.37|179.64|
|Mean net longwave W/m²|−104.94|−92.95|
|Mean net radiation W/m²|171.43|86.69|
|Source snow-positive hours|0|96|

Integrated second/first-window dose ratios are **0.72215 for positive degree days, 0.54485 for net shortwave and 0.42389 for net radiation**. This is a real forcing separation: an unfitted temperature-proportional model and a radiation-proportional model cannot both predict the same seasonal transfer. It is not yet evidence about which predicts ice ablation better. Reflected-longwave-corrected radiometric surface temperature can be computed assuming emissivity0.93, but is not an independent thermometer. Neither wind nor low RH proves actual evaporation; this archive lacks direct latent-flux or debris-water measurements, and precipitation is reconstructed.

## Conservation calculation: what would a missing process need to explain?

Let M be ice ablation in m water equivalent. The corresponding latent energy per horizontal area is E=ρw Lf M. Thus an interval residual ΔM implies a mean energy difference ρw Lf ΔM/Δt, if the footprints match. A generic **0.10m w.e. difference over30days requires12.886W/m²**. This is dimensional accounting, not an observed residual.

For a debris control volume, integrated ice heat is incoming net surface radiation plus turbulent/rain/advection energy, minus stored-energy change and other boundary losses. One melt residual constrains only this sum. It cannot identify its separate terms. Measured negative longwave cannot simply be subtracted from a trained PDD melt prediction: PDD calibration already absorbs some of those processes, and debris conduction filters their timing.

Using Lv=2.476MJ/kg, the generic12.886W/m² corresponds to0.450mm/day evaporation if every joule otherwise reached ice. This is a deliberately favourable upper-effect calculation; conduction and storage generally reduce instantaneous transmission. It demonstrates physical plausibility of an energy scale, **not a causal diagnosis**. A process explanation requires measured or independently bounded moisture and turbulent exchange, plus a consistently coupled heat budget.

## A bounded next question that could produce useful new evidence

**Does an independently parameterized thermal-memory model retain predictive information across a measured seasonal radiation–temperature decoupling, where a calibrated one-driver model does not?**

Freeze before scoring the later Pirámide interval: train-only parameter selection; source-derived h/k/porosity provenance and uncertainty; initial thermal profile/spinup; exact ice target; visit-time envelope; snow-treatment rule; source/instrument error semantics; professional-output component gate. Test a fixed factorial set: PDD, constant rate, net-SW, measured-net-radiation with debris storage, and full surface balance with the same conduction/geometry. A measured-radiation driver and a full balance have different information requirements and must be labelled separately. The input radiation is not a predicted flux.

Required primary result: later-interval native ice-ablation error, simultaneously over all declared clock/geometry scenarios, against identical target support. Supporting independent result: surface/debris temperature timing, if a genuine separate measurement exists. Ablate storage at matched training fit to test whether memory, rather than extra calibration freedom, causes improvement. Never select h, density, wetness, snow masking or clock by the later melt observation. Report a failure if every physically supported parameter set misses.

**Originality gate remains open.** This question is physically meaningful but its thermal-memory mechanism is established. A new transferable uncertainty certificate coupling energy, moisture stock and sensor/footprint support may be an engineering contribution only after an explicit prior-art audit and independent coverage/skill test. I do not call the conservation identity a new theorem or the proposed test a discovery.

## Fail-closed gates

1. Unknown professional CSV snow/ice variable → contextual output only, no certified professional-superiority comparison.
2. Snow-positive target period with no resolved snow treatment → conditional score; do not silently remove it.
3. Unknown initial thermal state and insufficient spinup → bound it or abstain; no fitted hidden heat reservoir.
4. Moisture/evaporation unobserved → a hypothesis or an independent-bound sensitivity, not a demonstrated cause.
5. AWS is200m from the stake and k is literature-derived → forcing/property transfer uncertainty remains; not co-located flux validation.
6. A single later interval and spatially related stakes → an event-level retrospective test, not consistent multi-glacier field superiority.

## Reproduce and verify

`python3 process_diagnostics.py ../targets/6_PIR_20142015.zip PIRAMIDE_PROCESS_DIAGNOSTICS.json`

`python3 -m unittest test_process_diagnostics.py test_resolve_target.py test_arolla_baseline.py`

All12 tests pass. They test dimensional energy accounting, signed flux, reflected downwelling radiation, incomplete support and invalid forcing, plus the earlier source/mapping gates. They do not validate environmental melt accuracy.
