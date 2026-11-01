# Independent public thermal-incision validation data audit

Bounded search completed 1 November 2026. No additional model was fit or scored. **No new fully matched, directly measured open-channel incision + water-temperature + hydraulic time-series CSV package was verified as downloadable in this search.** This is a bounded search outcome, not proof that none exists.

| Candidate | Public artifact actually verified | Relevant observations | Gate for independent test |
|---|---|---|---|
| Vincent et al.2010, Lac de Rochemelon, France | Publisher full HTML and PDF | Measured lake temperature, channel-bed elevations and geometry, plus lake-level-derived discharge during August2005 drainage | Strongest distinct surface-channel field event, but numerical source files not found; published graphs would require audited digitization and uncertainty |
| Pohle et al.2022, Rhonegletscher artificial moulins | Publisher full article and supplementary PDF; paper cites dataset DOI | Co-located pressure, temperature and tracer-derived Q/velocity/cross-sectional evolution | Different pressurised vertical-channel physics; exact archive file list/download not verified because DOI destination returned HTTP429 |
| Karlstrom et al.2014, Llewellyn Glacier | Publisher PDF and article download listing | Water-temperature/discharge monitoring, channel geometry and approximate surrounding-ice ablation | No paired directly measured channel-bed incision time series established; no separate numerical supplement listed |
| Isenko and Mavlyudov2002, Fisht Glacier | Society-hosted primary full PDF | Three four-day average measured incision values with channel widths and approximate slopes | No paired dated water-temperature/Q time series in the inspected field comparison; unsuitable for claimed exact hindcast |
| Isenko et al.2005 | Publisher abstract/full-paper landing record | Reports temperature observations from glacier streams | No downloadable matched numerical data package verified; abstract cannot establish paired timing/targets |

## Strongest potential distinct surface-channel test

Vincent et al., *Outburst flood hazard for glacier-dammed Lac de Rochemelon, France*, Journal of Glaciology56(195),91–100 (2010).

[Primary full article](https://www.cambridge.org/core/journals/journal-of-glaciology/article/outburst-flood-hazard-for-glacierdammed-lac-de-rochemelon-france/C183D075F0E9ABCD7CE5251D6096FCA4)

[Direct publisher PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C183D075F0E9ABCD7CE5251D6096FCA4/S0022143000206357a.pdf/outburst-flood-hazard-for-glacier-dammed-lac-de-rochemelon-france.pdf)

The paper documents 18 direct channel-floor surveys at cross-sectionsA–C, continuous lake thermistors with ±0.1°C accuracy, additional profile sectionsD–G, and a drainage start of25August2005 at1900h. Figures5/6 contain observed profiles and temperature/bed/discharge series; Figure8 combines model and measurements and must not be indiscriminately treated as observations. Original tabular time series were not located. Published curves could support an explicitly digitized conditional benchmark after uncertainty and spatial matching checks, but are not raw data. Cold ice is reported (−2°C at12m), and discharge is derived from lake level with recharge correction. The article explicitly states its heat-transfer coefficients were tuned to these observations, so using those fitted coefficients as an independent benchmark on this event would leak the target. This is a promising untouched event relative to current Swiss-model development, not a completed validation dataset.

## Distinct pressurised-channel lead: available citation, download unverified

Pohle et al., *Characterising englacial R-channels using artificial moulins* (2022), DOI10.1017/jog.2022.4.

[Primary article](https://www.cambridge.org/core/journals/journal-of-glaciology/article/characterising-englacial-rchannels-using-artificial-moulins/B47ADF6DAEC8D4B75497ABA0FABD1373)

[Public supplementary PDF](https://static.cambridge.org/content/id/urn%3Acambridge.org%3Aid%3Aarticle%3AS0022143022000041/resource/name/S0022143022000041sup001.pdf)

[Paper-cited data/code DOI](https://doi.org/10.3929/ethz-b-000504661)

The source explicitly identifies a raw-data/code deposit. A single redirected public request returned HTTP429 at the Research Collection destination; no retries or bypass were attempted. An exact archive filename, license, numerical coverage and bytes remain unverified. Tracer-derived cross-sectional area is not a direct bed-incision survey. This can potentially test englacial energy/area evolution after a separate applicability audit, not validate the present open cliff/stream system without adaptation.

## Other primary artifacts inspected

- [Karlstrom et al.2014 publisher record](https://tc.copernicus.org/articles/8/537/2014/) and [direct PDF](https://tc.copernicus.org/articles/8/537/2014/tc-8-537-2014.pdf): the record lists articlePDF, metadataXML and citation downloads, not a numerical data supplement. The field methods estimate bulk surrounding-ice lowering only at order-of-magnitude accuracy. Temperature/discharge observations should not be relabelled as independently measured channel incision.
- [Isenko and Mavlyudov2002 primary society PDF](https://web.seppyo.org/bgr/pdf/19/BGR19P93.PDF): the inspected Fisht comparison gives four-day average3.5/4.3/2.2cm/day incision with widths0.4/0.4/0.2m and approximate12° slope. No exact observation dates and paired water-temperature/discharge records are established. These are coarse physical plausibility checks, not a new high-accuracy temporal test.
- [Isenko et al.2005 publisher record](https://www.sciencedirect.com/science/article/abs/pii/S0165232X04001594), DOI10.1016/j.coldregions.2004.12.003: no independent machine-readable matched archive was verified.
- [Alexander et al.2021 primary preprint](https://eartharxiv.org/repository/object/2279/download/4689/): channel melt rates are calculated from velocity and an assumed0.01°C temperature difference, rather than directly measured melt targets. This fails the independent thermal-melt-target requirement even if drifter data are available.
- [SHED public model input/software](https://zenodo.org/records/7463322) and [model-region data](https://zenodo.org/records/7652297) are explicit accessible artifacts, but model outputs/forcing are not independent measured thermal-incision truth. They were excluded.

## Decision

Do not label any of these leads as newly attained consistent field superiority. The strongest bounded next public-data task is a provenance-preserving Rochemelon2005 graph extraction, with external-input and target curves separated, a model fixed before target digitization, all source uncertainty carried forward, and no reuse of Rochemelon-fitted heat-transfer coefficients as independent predictions. If the requirement is raw paired numerical CSVs rather than digitized published observations, it remains unmet in this bounded audit. No authors were contacted and no private data were used.
