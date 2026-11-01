# Reproduce the bounded Swiss thermal-component experiments

All files are additive. No active Claude source, frozen simulation or checkpoint is edited by these commands. Make a copy before running because each output directory is overwritten. Python3.9+ suffices for the solvers/audits/tests; figures additionally need NumPy and Matplotlib.

The original public source archive members must be extracted with their SOURCE_MANIFEST.json intact. Pass that directory as SOURCE below. Keep sibling temperature_digitization/ next to conditional_benchmark/. Source identity checks and original file hashes are preserved. Source temperature timestamps are naive while direct manual files state UTC; clock scenarios retain that unresolved issue.

```sh
python -m unittest discover -s conditional_benchmark -p 'test_*.py'
python -m unittest discover -s temperature_digitization -p 'test_*.py'
```

There are4 closure checks and15 temperature/gap checks, all passing in the final portable bundle.

From conditional_benchmark/:

```sh
python run_comparison_v3_executed.py SOURCE --out NEW_STRICT_RESULTS
python run_comparison_v3_executed.py SOURCE --forcing author_daily --out NEW_DAILY_RESULTS
python run_v4_native_gaps.py SOURCE --out NEW_V4_RESULTS
python run_v5_local_geometry.py SOURCE --out NEW_V5_RESULTS
python audit_comparison.py NEW_DAILY_RESULTS/CONDITIONAL_RESULTS.json
python audit_comparison.py NEW_V4_RESULTS/linear/CONDITIONAL_RESULTS.json
python audit_v5.py NEW_V5_RESULTS NEW_V4_RESULTS
python make_figures.py NEW_DAILY_RESULTS/CONDITIONAL_RESULTS.json --out NEW_DAILY_FIGURES
python make_figures_v4.py NEW_V4_RESULTS/linear/CONDITIONAL_RESULTS.json --out NEW_V4_FIGURES
python make_figures_v5.py NEW_V5_RESULTS/linear/CONDITIONAL_RESULTS.json --out NEW_V5_FIGURES
```

For V4 arithmetic audit repeat for left and right output folders. The strict run intentionally produces no field scores because all complete manual intervals have unsupported native temperature. It must not be passed to an audit expecting scored scenarios. run_comparison_v3_executed.py exactly matches the source hash in the saved author-daily results; run_comparison.py adds the explicit temperature-factory API used by V4. All underlying numerical model families remain unchanged.

V4 gap methods are fixed from a separate temperature-only mask experiment. V5 changes only the imposed local-geometry construction, without station-specific learned gains. V5 was selected after earlier target errors were exposed and is labelled development. Reports preserve failures against trained constant-rate/h and unfitted heat-transfer comparators. Engineering scenarios are not confidence limits or independent replicates. No complete coupled cliff-model field accuracy is asserted.

Read RESULTS.md, RESULTS_V4.md, RESULTS_V5_DEVELOPMENT.md and HYDRAULIC_INPUT_QUALITY.md before using the results in a paper. The V5 common-cohort audit controls both training and test membership when isolating the geometry feature. Source daily hydraulic depths remain pressure-derived proxies rather than independently established continuous local depth observations.
