# Matched daily melt scoring

`paired_skill.score` requires the frozen full site/date support and complete predictions from every declared comparator. It refuses nonfinite data, duplicate records, missing predictions and mixed units. It reports pooled and equal-site MAE, signed bias, RMSE and descriptive day-win counts. Seven hand-calculated/refusal tests pass.

No favorable subset is selected automatically. Replacing the expected support requires a versioned protocol, not silent omission of difficult cases. Common units alone do not establish matching observation support: interval lengths, snow/refreezing definitions, sensor footprints, prediction timing and training independence must be verified by the data contract before scoring. The function does not certify any of them, compute a confidence interval or establish historical novelty.
