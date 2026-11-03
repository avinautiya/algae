# Conditional plate-signal experiment — prospective scoring protocol

## Claim and native-data gate

This is a conditional diagnostic of archived signed raw plate readings under observed boundary temperatures. It is NOT prediction of melt, latent heat, cliff retreat, future weather, material properties, or a new thermodynamic law. The archived hourly sensor records are means; their exact averaging interval and native plate orientation/position are not independently resolved. Therefore the **native interval-matched flux-validation gate is NO-GO**. If executed, the separate diagnostic explicitly treats each recorded mean as a representative point at its labelled timestamp. Scores cannot be called exact interval-flux validation. This convention must not be changed after seeing results.

Primary data: Changri Nup archive DOI10.5281/zenodo.3048780, verifiedMD5 4797deeda5c93a1895f08efb42183086. Input is the source-aligned hourly audit from the verified loader, not SR50 ablation. Measurements: Hukseflux interface plates at nominal10cm/site1 and8cm/site2; actual sensor positions may migrate. Top probes were exposed/replaced on1May2016; this protocol stays before that change. Frozen temperatures remain allowed because the mathematical problem is a sensible subcolumn with MEASURED bottom boundary; no conversion of flux to ice melt is authorized. A constant effective-conduction hypothesis may fail due phase change, advection or sensor placement.

## Dates, fitting and baselines

Initialize once at the FIRST archived source row (28Nov2015). Train:1Dec2015inclusive–29Dec2015exclusive, Nepal UTC+05:45. Test:29Dec2015inclusive–12Jan2016exclusive. Do not reset at either boundary. The first14source days were already used for training-only diagnostics; that exposure is disclosed. January targets have not been scored when this specification is frozen.

Site1 uses measured1/5/10cm probes. Fit alpha on training interior5cm temperature only,61log-grid values1e-8–1e-5m²/s. The1cm and10cm temperatures are imposed observed boundaries. Site2 uses1/8cm observed boundaries and inherits site1alpha because it has no interior probe; its initial profile is linear between first observed endpoints.

For each site predict basal gradient from the64mode Dirichlet solution. Fit one through-origin raw plate coefficient from TRAINING ONLY. The coefficient combines effectiveconductivity magnitude and inferred orientation; neither is a separately validated intrinsic property. Negative coefficient is retained as training-inferred orientation, never flipped on the test set. No published full-season optimized conductivity enters fitting.

Comparators: each site's training median; instantaneous measured-endpoint gradient with its OWN training coefficient; causal observed-endpoint ridge using lags0/1/3/6/12/24hours. Ridge penalty0/.01/1/100 selected by the final7training days after fitting on first21days, then refit on all28training days. The runner refuses non-hourly or gapped endpoint support before building hour-labelled lags. Ridge receives only the observed endpoint columns. Each comparator predicts exactly the same finite-plate test timestamps. Missing plate records remain missing; no favorable target subset selection.

Primary contrast: diffusion versus causal ridge, equal-site mean MAE. Pass only if diffusion is better at BOTH sites, aggregate MAE is at least5% lower, both numerical gates pass, and fitted alpha is not a grid endpoint. This rule does not establish statistical significance. Median and instantaneous-gradient comparisons are secondary. The opposite fitted plate-orientation sign is reported as a sensitivity, never chosen from held-out errors.

Fit alpha once for site1, not each testsite. All61training objectives and both sites' parameters must be retained. Alpha endpoint selection is reported. No superiority promise and no model winner selection affect the protocol.

## Numerical and scoring rules

Report every site's n,MAE,RMSE,signedbias in rawW/m². Day/hour records are serially correlated; do not label them independent field replications. No confidence interval or climatic transfer claim follows.

Check64vs128mode instantaneousgradient predictions with the SAME fitted coefficients. Numerical gate: max heldout signal difference must not exceed max(0.5W/m²,1%trainingplateRMS). Retain scores even when this gate fails, clearly invalidating accuracy claims. A native averaging-aware test needs new documented alignment and a new freeze.

The code writes a fresh output directory only. Verify all frozen source/data hashes and publish the freeze commit BEFORE scoring. A later correction must preserve the original source/results and record the amendment; it must not silently rewrite this protocol.

## Prior art and novelty ceiling

The diffusion/conduction inversion and flux-scale identification are established inverse heat-transfer ideas. DCG-MIP2026 identifies conductivity estimates, parameter provenance, model process attribution and colocated observations as unresolved practical issues (https://doi.org/10.5194/tc-20-1895-2026). Collier2014 and Giese2020 already model debris moisture/phase reservoirs (https://doi.org/10.5194/tc-8-1429-2014 ; https://doi.org/10.5194/tc-14-1555-2020). Neither moisture nor conduction is claimed new. This experiment could test whether independently trained plate-constrained parameters transfer better than a direct baseline under a declared observation approximation. It does not establish a novel glacier-melt predictor or professional-model superiority.
