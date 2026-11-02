# Frozen field-temperature experiment

Read `REPORT.md`, then `results/FIELD_TEMPERATURE_RESULTS.json`. Inputs and protocol are in adjacent `source/` and `protocol/` folders. Saved prediction CSVs include raw observation values and role labels; individual profile JSONs include all alpha-grid objectives and fitted lag parameters. `RIDGE_PARAMETERS.json` stores the training-only coefficients.

Reproduction (Python with NumPy/SciPy/Matplotlib):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/python3 run_benchmark.py --out results_reproduction
/usr/bin/python3 audit_saved_results.py
/usr/bin/python3 -m unittest test_leakage.py
/usr/bin/python3 make_figures.py
```

The original result directory cannot be overwritten by the runner. For byte-identical source provenance use a separate copy of the package with `solver/dirichlet_slab_temperature_v1_executed.py` copied to `solver/dirichlet_slab.py`; the current solver has an additive gradient API. Audit and plotting scripts read the delivered original `results/` directory. Runtime and platform metadata may differ on reproduction.
