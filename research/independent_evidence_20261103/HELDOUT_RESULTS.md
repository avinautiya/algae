# Frozen conditional test: no superiority over the strong baseline

Source, protocol and data hashes were published BEFORE scoring at commit `2a0828dd2a4c8c42562a277eec8dfb2b991cf035`, prefix `research/independent_evidence_20261103/`. The conditional representative-point diagnostic then ran once, using the predeclared28-day training/14-day evaluation windows. All672 fixed evaluation records are retained. There was no parameter change after scoring.

| Model | Equal-site MAE, raw W/m² |
|---|---:|
| Causal boundary ridge |0.8465|
| Diffusion with training-only flux scale |1.1664|
| Instantaneous gradient with training-only flux scale |3.8885|
| Training median |9.5479|

The diffusion model is substantially better than instantaneous-gradient and constant predictions, but worse than ridge at BOTH sites. The predeclared superiority rule therefore FAILS. This is not rescued by numerical error:64-versus128-mode signal differences are0.0231 and0.0185W/m², below the frozen0.5W/m² threshold; the fitted diffusivity is interior to the grid at3.5481e-7m²/s. Those numerical checks do not prove the physical model is adequate.

The result is conditional on measured boundary temperatures and a representative-point reading of archived hourly means. Native averaging intervals, orientation and sensor position are not independently resolved. Site2 inherits diffusivity but fits its own amplitude locally. The two sensors belong to one glacier. This is neither zero-shot geographical transfer nor a weather forecast, and no latent-melt, full-cliff or climate-model accuracy follows.

The primary outcome is negative against the strong comparator. Retain the result, plots, complete predictions, training searches and sign sensitivity. Do not reuse this inspected window as a fresh validation set for a modified model. Any subsequent candidate needs a newly justified mechanism and a separately declared evaluation design.

An independent post-score audit verified the remote pre-execution freeze and source blobs, matched every raw target and the336+336 fixed support, recomputed all metrics to3.55e-15 and replayed fixed-parameter predictions to7.55e-15 without refitting. Diffusion MAE is37.79% higher than ridge. See `decision/POSTSCORE_AUDIT.md` and its executable audit.

The professional baseline calculation is a different comparison:−21.975864% from a conditional raw-calendar reconstruction versus the published−23.3%. Its unresolved1.324136percentage-point discrepancy is NOT an improvement by a new model.
