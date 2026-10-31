# Second-pass source and data audit

Checked 31 October 2026. This supplements the previous packet's SOURCES.md. No exhaustive-priority guarantee is made. Focused searches covered open-channel thermal erosion, hydraulic-jump heat exchange, nonlinear ice morphology, CFD/phase-change implementation and predictive validation data.

| Source | What was inspected | Consequence |
|---|---|---|
| [Jarosch & Gudmundsson 2012](https://doi.org/10.5194/tc-6-493-2012) | Full open PDF downloaded; §§2–3 and discussion inspected, equation page visually checked. | Channel evolution is established; water hydraulics use Manning and local melt distribution is prescribed. Explicit water/heat resolution was identified as a possible extension. |
| [Perissutti et al. 2025, JFM](https://doi.org/10.1017/jfm.2025.10615) | Full publisher methods, scaling, cost and conclusion passages. | DNS/phase-field ice morphodynamics and scale extrapolation already exist. Geometry is ice above a stream, not a free-surface supraglacial step train. Physical melting timescale matters; the paper uses inverse-Stefan convention. |
| [Perissutti et al. 2024](https://doi.org/10.1016/j.ijmultiphaseflow.2024.105007) | Publisher abstract/results summary. | Existing turbulent ice-pattern simulations rule out a generic first-CFD or first-pattern claim. |
| [Bushuk et al. 2019 open manuscript](https://www.gfdl.noaa.gov/wp-content/uploads/2019/05/Bushuk_JFM2019.pdf) | Primary abstract/conclusion passages. | Ice geometry–turbulence–heat feedback and migrating equilibrium scallops are established. Scallops are not the same as free-surface jump-bounded steps. |
| [Ogier et al. 2021](https://doi.org/10.5194/tc-15-5133-2021) | Full methods, thermodynamics, measurement uncertainty and data-availability passages. | Potential field evaluation of incision and cooling; no proof that a heat-transfer discrepancy is caused by steps. |
| [Ogier data collection](https://www.research-collection.ethz.ch/entities/researchdata/27f98ac0-2289-4539-9fb6-c6ee6308600a) | Public archive listing and licence. | MIT-licensed record lists 322.20 MB ZIP and 847 B README. Raw contents not inspected: downloads returned 403/429/cache errors. |
| [Mukha et al. 2022](https://doi.org/10.3390/fluids7030101) | Primary methods/results summary. | Hydraulic-jump LES is sensitive to interface capture, mesh and numerical dissipation. Hydraulic jump benchmarks do not validate melt flux. |
| [Modelling the evolution of elliptical Röthlisberger channels, 2026](https://www.cambridge.org/core/journals/journal-of-glaciology/article/modelling-the-evolution-of-elliptical-rothlisberger-channels/981B8BCE9C395330DDC3736ABB045F69) | Publisher abstract screening; full method not audited. | Evolving channel shape and nonuniform melt are active prior art; this concerns enclosed channels, not our initial open reach. |
| [Karlstrom et al. 2013](https://doi.org/10.1002/jgrf.20135) | Primary abstract/method screening. | Curvature-dependent melt and supraglacial meander instability already exist. Do not add planform meandering as an unexplored novelty claim. |
| [OpenFOAM standard solvers](https://www.openfoam.com/documentation/user-guide/a-reference/a.1-standard-solvers) and [melting source](https://api.openfoam.com/2606/classFoam_1_1fv_1_1solidificationMeltingSource.html) | Official capability documentation. | Useful components, not proof of turnkey coupled ice/free-surface capability. |
| [Basilisk embedded boundaries](https://basilisk.dalembert.upmc.fr/src/embed.h) and [solidification/melt method](https://www.sciencedirect.com/science/article/pii/S0021999122008920) | Official module and primary method screening. | Alternative numerical building blocks; coupling and verification remain work. |

## Observational constraints that change the design

Ogier's measurements are from one Swiss glacier drainage event, not independent Himalayan validation. The paper reports 11 longitudinal surveys with roughly 0.1–0.5 m vertical uncertainty. Temperatures are point measurements near the channel bottom, with nominal accuracy ±0.05 or ±0.1 °C depending on station. Bulk temperature can differ systematically. Daily wall marks are interpreted as incision proxies. These distinctions belong in an observation operator and uncertainty analysis; do not directly score a bulk-temperature prediction against a near-bed sensor or treat every logger timestamp as an independent test.

The field paper's discrepancy with generic heat-transfer correlations is motivation, not attribution to a particular mechanism. Sparse geometry, sensor placement and closure uncertainty can produce competing explanations.

## Raw-data readiness

Confirmed: record identity, archive/README listing, stated licence and article measurement descriptions.

Not confirmed: internal filenames, time-axis conventions, sensor calibration files, initial full geometry, recoverable spatial temperature fields, mesh-ready topography, quality flags or completeness of boundary histories. These require the archive. No raw measurements were used to fit or score a model in this session.

The paper on cyclic ice steps links original photos at `https://sites.google.com/site/miwayokokawawebsite/open-data`; this link could not be retrieved here. Printed figures/tables remain available, but digitization has uncertainty and is not equivalent to recovering original image sequences.

## Sources that still require full priority checks

Finish forward-citation searches for the closest Wu/Izumi, Yokokawa and Camporeale papers; inspect current open-channel CFD/thermal-erosion work under alternative terminology. Read the full 2026 enclosed-channel paper before making any cross-sectional-shape claim. Maintain the previous exclusions: no first thermal ice model, first memory model, first moving-boundary CFD, first spectral method or first ice-step discovery.

One new opening may be the particular combination of resolved free-surface jump heat transfer, finite-amplitude ice-step stability and transferable finite-reach predictions. It remains a candidate combination until its mechanism and benefit are demonstrated; combining established ingredients is not automatically a major discovery.
