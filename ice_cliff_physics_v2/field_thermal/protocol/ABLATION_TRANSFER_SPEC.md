# Conditional ablation transfer: separate from temperature-only validation

Written before ablation model scoring. Source mass totals and metadata are already inspected, so this is retrospective spatial transfer, not blind discovery. Do not use published profile conductivities as independent predictors: those estimates may have used the same ablation data.

## Explicit author-linked targets only

Four profile/stake links are admissible: KDT_2020b/P_AWS_1 and KDT_2023a/P11,b/P12,c/P13. The association is documented by the source, not inferred from nearest coordinates. They are nearby linked locations, not proven centimeter-identical material columns. Temperature and stake debris thicknesses differ; P12 is particularly important (profile32cm versus stake29cm).

Use only the summer ice-ablation record for each link. Do not reuse the2020b following winter balance: its thermal forcing is outside the record. Source negative mass change means positive loss; retain the original kg/m2 unit. No conversion to water-equivalent depth or geometric ice volume is made without an explicit density.

Ablation dates lack hours. Primary clock scenario is00:00 source ADT on each date (UTC=ADT+8h), a declared convention rather than established survey time. Also report every Cartesian start/end shift in {-12,0,+12}hours. No clock is selected by target error.

## Distinct full-column model

The primary temperature test remains bounded between measured probes. Ablation needs a different physical boundary: measured shallowest temperature at ztop, assumed temperate ice at depth d with T=0°C, and all remaining real sensors within that domain used only in the early temperature-calibration window. This bottom is a phase assumption, not an observed thermometer. `Base`0 in2020b is still not an observation.

Run geometry scenarios independently with d equal to author profile debris thickness and linked stake thickness. A real sensor beyond a scenario's assumed interface cannot be silently used as a debris observation; exclude it from that scenario by geometry before fitting and retain the consistency flag. Do not tune depths or choose which thickness gives a favorable loss prediction. Fixed depth and a common material within each column remain explicit assumptions.

For2020b the original raw sensor temperatures cover the full summer ablation dates, but installation configurations differ between the older June and later July setups. The transition date is not established. Use confirmed processed geometry from27July onward; the first~2days of the mass interval require explicitly separated earlier/current depth configurations (upper0.5cm versus5cm, with corresponding second sensor18cm versus15cm). This resolves missing temperature, not missing geometry or survey hours.

For2023 the temperature series starts10June20UTC after the8June ablation start. The latest archive subset has the same thermal records, so no earlier data have been verified. Do not label the full season100% observed forcing.

## Model families and calibration

Use the same physically justified boundary inputs for all families:

1. Distributed normal diffusion with effective alpha fitted from first14days of temperature only, with the same predeclared bounds/initialization rules as appropriate to the full-column geometry.
2. Quasi-steady conductive profile driven by the same shallow temperature and basal0°C boundary.
3. A one-node normal diffusion approximation, with top and bottom resistances2k/L and node capacity Cv*L; its relaxation is L²/(4alpha). Fit its effective alpha from the same early temperature objective, not ablation tests.
4. Training-origin constant daily ablation, using the source stake's full positive loss and documented interval duration. This is a weaker statistical baseline and cannot be the sole comparator.

For each thermal family, predictions before conductivity calibration are a unit-conductivity integrated basal-gradient quantity J. Only this diffusivity/shape response is determined by temperatures. Fit one nonnegative effective conductivity separately for each family on the SAME source ablation target. Do not give a baseline the conductivity optimized for the PDE and call that equal calibration. Keep its implied Cv=k/alpha as a derived effective value, not an independently measured property.

Primary transfer matrix: every one of the4source stakes calibrates its own scalar; predict each of the other3 targets, giving12source-to-target comparisons. Never select a successful origin and suppress the others. These are not12independent glaciers. A separately labelled four-fold leave-one-stake-out fit on all3other stakes may be added, provided all folds and families are reported; do not conflate that with single-origin transfer.

## Missing period and initial state

Keep the measured forcing support and predicted loss over that support separate from the full stake interval. Initialize the covered-period temperature state from measured real probes at the first common thermal timestamp, plus the declared basal0 boundary. No backward diffusion or inferred unobserved interior history is permitted.

For a SECONDARY full-period point estimate, fill only the missing leading loss by periodically repeating the model's first complete observed24h unit-conductivity basal-flux pattern, matched by source time of day. This is a model-based flux reconstruction, not newly measured forcing. The full-period initial-state uncertainty is not solved by this procedure.

Report leading-loss sensitivity as0and as the maximum24h unit-conductivity positive basal heat among the first14observed days times missing duration/24h. The latter is an empirical engineering upper scenario, not a proved meteorological bound. No later target loss or winning model error chooses the pattern. If a whole observed day is unavailable, do not invent one.

If source conductivity is calibrated with an uncertain missing leading contribution, retain its corresponding interval: for positive unit-gradient integrals Jcovered and Jmissing,max, k lies between Lf*Msource/(Jcovered+Jmissing,max) and Lf*Msource/Jcovered under these assumed nonnegative-ablation scenarios. This is not a statistical confidence interval. When original source forcing/geometry is complete, no artificial gap term is added.

## Primary necessary check versus secondary point accuracy

The covered-period predicted positive loss must not exceed full measured summer loss if ablation is nonnegative in the unobserved period and there is no omitted accumulation/refreezing. This is a ONE-SIDED consistency/falsification check under those assumptions. Being below the full measured loss does NOT establish covered-period accuracy, because no covered-only mass target was measured.

Full-season absolute/relative errors are SECONDARY, explicitly conditional on the leading reconstruction, date clock and thickness scenario. Score every fold, model and scenario; no favorable scenario is selected. Retain negative basal heat separately as cooling/refreezing demand. A melt-only positive-flux approximation must be labelled; do not silently convert signed heat into measured net mass balance or claim phase physics is solved.

## Numerical and scientific checks

Verify flux/phase-energy convergence independently of temperature convergence. For homogeneous constant material and two Dirichlet boundaries, signed basal heat obeys

    Qbottom = integral[k(Ttop-Tbottom)/L]dt - Cv*Delta integral[(z/L)T]dz.

Thus long-period diffusion versus quasi-steady differences may mostly be endpoint storage, even when their temperature dynamics differ substantially. Do not promise seasonal PDE superiority from better interior temperatures. Positive-flux clipping is a separate nonlinear approximation and must not be hidden in this signed identity.

A successful result supports conditional sub-debris ablation transfer with local thermal calibration, not a field validation of water-film cooling, stream routing, cliff collapse or Himalayan glacier-wide prediction. Report the actual input/observation mismatch and uncertainty rather than treating these point targets as full-cliff truth.
