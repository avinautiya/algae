# Development decisions and failed cases

- The first downstream-wetting connection used upstream pool stage. Independent
  review rejected that identification. The current model requires explicit
  downstream Manning roughness and friction slope and reports both stages.
- A necessary head-admissibility check was added. It does not solve backwater.
  The synthetic example weir coefficient changed from 1.7 to 1.2 m^0.5/s to
  retain head for the explicitly modeled contact reach. This was a development
  admissibility change, not calibration to observations or a held-out outcome.
- The high-friction-slope variant remains inadmissible and is reported as
  refused, not removed from the development study. All six attempted variants
  and their configurations are retained in `development_results.json`.
- Extremely short forcing intervals initially could be skipped by a 1e-10-s
  stopping tolerance. Strict advancement and a regression case replaced that
  behavior. No production campaign used this development implementation.
- Surface reservoirs use extensive states with explicit new-area material.
  Separate volumetric cold cost accounts for newly melted cold ice; this is an
  effective enthalpy model, not a resolved cold-ice conduction PDE.
- Input weather remains observed station data with flags. No automatic
  station-to-face radiation transformation, RH clipping, missing-rain filling
  or air-to-water temperature substitution was performed.
