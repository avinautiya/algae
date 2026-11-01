# Environmental extension — evidence and integration record

This additive package leaves Claude's active worktrees, frozen campaigns and
checkpoints unchanged. It supplies new code and genuine Himalayan geometry;
it does not retroactively turn the earlier plan-view runs into cliff forecasts.

## What the new environmental question is

Does a model that predicts the timing and position of stream–ice contact give a
more useful forecast of **exposed cliff area and melt over a common 3-D footprint**
than a model that prescribes contact from debris area or removal time?

Test this over nonzero under-rock/lateral transmission, environmental heat
forcing and uncertain geometry. A result is consequential only if its effect
exceeds numerical and closure uncertainty and improves the stated prediction
against an independent reference. If transmission or atmospheric forcing erases
the mechanism, report that rather than designing a range that preserves it.

This is a candidate contribution. Adding established energy balance, polygon
geometry or a porous rating is not itself a new mathematical or environmental
discovery. Existing cliff models already include several of these processes.

## Genuine Himalayan data now acquired

Kraaijenbrink and Immerzeel's [Zenodo release 1.2.0](https://zenodo.org/records/13073856)
contains Lirung and Langtang cliff-profile point clouds. Selected members were
retrieved with HTTP byte-range and size checks; local SHA256 values are recorded
in `data/profile_download_manifest.json`. The full 463.9 MB archive MD5 was not
verified because only selected members were downloaded. The source is CC-BY4.0.

The CSV actually contains **391,911 points, 11 profiles, 83 profile–survey
combinations and 9 cliff identifiers**. All profiles are included in the
reproducible figure. These are corrected, unordered terrain point clouds, not
classified cliff/stream/debris masks. They do not supply coincident stream
discharge, temperature or collapse labels. Binning does not remove survey error.
The accompanying survey table gives non-reference vertical RMSE values from
0.1817 to 0.6194 m. Its October 2015 reference rows contain zeros by construction;
they are not evidence of error-free terrain. Survey-wide relative RMSE is not
a cliff-normal uncertainty or an independent confidence interval. In particular,
these profiles cannot independently resolve millimetre-scale contact thresholds.

The [NASA Langtang AWS product](https://doi.org/10.5067/TNQV3XTMP9CP) provides
meteorology for 2017–2018. Radiation, stream runoff and stream water temperature
must be checked separately; air temperature must not be relabelled as water
temperature. The forcing agent's source audit records file-level accessibility.

## Physical and evidential distinctions that must survive integration

1. **Floor incision is not cliff-normal retreat.** `integration/stepped_cliff.py`
   represents both directions in an explicit extruded vertical-step geometry.
   Corner removal is counted once. The whole face retreats uniformly, so this
   geometry cannot generate undercuts or fracture. Its fixed crest and planar
   surfaces are assumptions.
2. **Heat is assigned once.** Advected heat must be depleted as it warms cold
   ice, enters debris storage, melts the floor and reaches the face. Face
   atmospheric energy cannot be counted as floor incision. Geometric ice
   removal from collapse cannot be counted as latent melt or runoff.
3. **Transmission must preserve the actual ice floor.** A below-debris gap may
   admit water below the debris crest, but not below an intact floor. A lateral
   or englacial route to a lower cliff elevation is a different geometry and
   cannot be represented by an arbitrary nonzero leak coefficient.
4. **Hydrostatic diagnostics are rejection criteria, not validation.** Froude
   number alone does not justify a hydrostatic closure at a steep crest. Curvature,
   aspect ratio, flow change and boundary conditions need independent checks.
5. **Cold ice and debris storage introduce energy states.** A negative cold
   reserve has to be paid before latent melt; warm debris may later release
   heat. Algebraic instantaneous resistance is not a storage model.
6. **Collapse cannot be closed by a guessed threshold.** A prescribed collapse
   perturbation can test sensitivity, but predicting collapse requires a
   defensible mechanical model and an independent verification target.

## When the manuscript may replace a limitation with a result

An implemented component can be described in Methods after its component
checks pass. A coupled environmental claim requires a compatible coupling,
water/heat/phase/mass ledgers, resolution and timestep convergence of the
prediction, and a comparison using the same ice footprint. A field prediction
also requires the actual site inputs. A sensitivity interval is not a confidence
interval, and synthetic forcing is not environmental validation.

Until those gates pass, the main manuscript's limitations must remain accurate.
The package is designed to replace specific missing components with tested
methods, not to delete uncertainty by wording.
