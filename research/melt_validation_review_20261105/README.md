# Glacier ablation transfer audit — pre-publication evidence release

This is a reproducible external-validation package, not a claim of a universally superior model or a completed journal submission. It preserves a new camera-stake test of previously frozen recipient-excluded predictions and published-model comparisons. Start with followup/FOLLOWUP_REPORT.md.

The candidate had2.865mm ice/day equal-site MAE on112complete stake-days at three locations in Pirámide's2022–2023campaign. It beat published KP1 fixed effects but did not beat every simple control under pooled-day weighting. On97days with common lagged-forcing coverage, it lost to all four DETI equation reconstructions. DETI variants are not certified author-code reproductions. The earlier Djankuat failure is retained in the accompanying report and must not be omitted from a paper.

## Reproduce

Tested versions and platform are in ENVIRONMENT.json. Use Python3.9or later compatible with the pinned requirements. No GPU is required. Installation time has not been measured on a clean machine.

```sh
python3 -m pip install -r requirements.txt
python3 download_source.py
python3 followup/piramide_test/audit.py
python3 reproduce.py
MPLCONFIGDIR=/tmp/glacier-release-mpl python3 followup/piramide_test/plot.py
```

download_source.py retrieves the62.7MB official CC BY4.0dataset and verifies its published MD5/size. reproduce.py runs scoring in a temporary directory and checks the saved score tables, protecting the released results from overwrite. The native scorer refuses existing result filenames. The source table decoder is a narrow validated reader of the archive's MATLAB5MCOS layout; no third-party MATLAB code is executed.

## Scope and provenance

The source is Melo-Velasco et al., Method Dependence in Thermal Conductivity and Aerodynamic Roughness Length Estimates on a Debris-Covered Glacier, https://doi.org/10.1029/2025JF008360 ; data https://doi.org/10.5281/zenodo.11581073 . Cite both and the original data creators when reusing derived observations. See DATA_LICENSE.md.

Published comparator references are Carenzo et al.2016, https://doi.org/10.1016/j.advwatres.2016.05.001 ; Winter-Billington et al.2020, https://doi.org/10.1017/jog.2020.57 ; and DCG-MIP2026supplement, https://doi.org/10.5194/tc-20-1895-2026-supplement . Their papers and software are not redistributed as our work.

Models contain no Pirámide calibration data. The local coefficient freeze preceded inspection of the new target bodies; endpoint eligibility was operationally resolved during the observation audit, so this is not externally registered blind confirmation. Source datetimes are unzoned calendars. Both hourly-label interpretations are retained. Earlier UTC-suffix serialization labels were corrected without changing any numeric score; historical copies remain in the local project, and the correction is documented in the report.

Training coefficients and method specifications are provided, but this archive does not redistribute the older donors' raw archives or derived training table because their redistribution terms were not explicitly available in the checked metadata. Full raw-source reconstruction of training remains separate work. Do not call this a complete training-data archive or DOI-backed repository deposit.

The author order supplied by Avi is Avi first and Dr. Manish Mehta (Wadia Institute of Himalayan Geology) coauthor. This pre-publication release does not assert final manuscript approval, contribution certification, funding or conflict declarations. Those require author confirmation.

OpenAI Codex and Anthropic Claude assisted code, methodological exploration, debugging and drafting during the project. Human authors must verify results, sources and interpretations before submission. AI systems are not authors. No journal acceptance, major historical novelty or Nature-level significance is asserted.
