# Arolla measured ice-ablation evaluation

This package advances from sensor-signal fitting to a comparison with actual stake-measured ICE ablation. It is a small retrospective temporal holdout, not a blind prospective trial, a new physical law, or validation of the stream/cliff solver. Source targets were inspected during data/units auditing before code was written. No evaluation interval fits a model coefficient or selects a timing scenario.

## Actual field target

The original Arolla archive DOI10.5281/zenodo.3047649 supplies cumulative stake readings and melt directly in mm water equivalent, using900kg/m³. H1D is the source-designated nearest stake to the AWS with known6cm debris depth. Its first interval supplies846mmw.e. for calibration; its second supplies315mmw.e. for evaluation. Three other stakes are secondary and share the same weather event, not independent glaciers. Dates are known but visit hours are not. The weather archive encloses all proposed visit hours.

The metadata says uncertainty0.1m, without a probability distribution or a declaration of whether it applies to readings or interval differences. It is not a Gaussian sigma or95%confidence interval. Both explicitly labelled bounded interpretations are propagated. The middle visit is shared between calibration and evaluation, so its measurement/timing error is not duplicated independently.

Snow is unobserved in the archive after4July; no code substitutes zero. Empirical models omit explicit snow processes, which can affect their errors. Conductivity is not used to fit them. A table-flattening error attributing Miage borrowing to conductivity was corrected before predictions; the note actually belongs to roughness. Canonical source provenance and the amendment are retained.

## Frozen comparison

`model/PROTOCOL.md` and `model/FREEZE.json` fix sources, data, models and reporting rules BEFORE predictions are computed. Three established one-coefficient local baselines use only the first interval: constant-rate persistence, positive-degree-day and net-shortwave-dose. The experiment does not relabel these as a novel model. All972 native-target model/timing rows remain visible. H1D is primary; other stakes are secondary and use the same AWS forcing as an approximation.

Actual published DEB_CF Arolla output is included for81 H1D contextual timing comparisons, without local rescaling. Its exported ice-only versus ice-plus-snow component is unresolved. Therefore these differences are NOT automatically like-for-like pure-ice scores and cannot establish superiority over professional models. Provided-property independence is also not fully documented. The source run, column, member checksums and clock assumptions are recorded.

Each method reports27 shared-visit-hour combinations under3 hourly-bin conventions:81 sensitivity scenarios, not81 independent observations. The nominal convention was declared before predictions. No scenario is selected by error. The timing helpers provide independent monotonicity and shared-reading error checks; neither is claimed historically new mathematics.

## Interpretation

This can establish quantitative retrospective agreement, or disagreement, with native ice-ablation measurements under stated assumptions. It cannot establish generalization across glaciers, a novel predictive mechanism, exceptional field accuracy, complete physical modelling, or Nature/competition outcomes. A favorable point error at one primary stake is insufficient for those claims.

The Changri Nup reproduction audit is retained separately:900kg/m³ was prescribed for submitted model outputs, not proven for author processing of every observed surface-height series. Density or measurement operators must not be selected merely to reproduce a desired score.
