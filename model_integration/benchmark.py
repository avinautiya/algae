"""Score precomputed, independently held-out predictions with cluster bootstrap.

This does not train models or certify independence from a manifest alone. The caller
must freeze complete calibration/preprocessing on training data before predicting.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .coupling import _array

REQUIRED_MODELS = ("no_algae", "established_algae", "empirical_satellite", "pigment_cell")


def audit_split(training_ids, test_ids, training_groups, test_groups):
    """Groups should be site-year or independent spatial/temporal blocks, not rows."""
    if not test_ids or not test_groups:
        raise ValueError("test samples and independent groups required")
    if len(set(test_ids)) != len(test_ids):
        raise ValueError("duplicate test sample IDs")
    if set(training_ids) & set(test_ids):
        raise ValueError("training/test sample leakage")
    if set(training_groups) & set(test_groups):
        raise ValueError("training/test group leakage")


def score(observed, predictions, groups, *, bootstrap=2000, seed=0):
    """Identical samples for all models; paired RMSE improvements, cluster resampling.

    RMSE/MAE/bias are observation-weighted. Bootstrap resamples whole groups;
    intervals describe sampling uncertainty, not calibration/model uncertainty.
    With <2 groups no interval is reported; even 2 groups is weak external evidence.
    """
    y = _array(observed, "observed")
    g = np.asarray(groups, dtype=str)
    if y.ndim != 1 or not len(y) or g.shape != y.shape or np.any(g == ""):
        raise ValueError("nonempty 1D observations and group IDs required")
    if set(predictions) != set(REQUIRED_MODELS):
        raise ValueError(f"must compare exactly {REQUIRED_MODELS}")
    p = {k: _array(v, k) for k, v in predictions.items()}
    if any(v.shape != y.shape for v in p.values()):
        raise ValueError("all models must predict the same observed samples")
    if bootstrap < 100:
        raise ValueError("use at least 100 bootstrap replicates")
    unique = np.unique(g)
    def metrics(idx):
        return {k: dict(rmse=float(np.sqrt(np.mean((v[idx] - y[idx]) ** 2))),
                        mae=float(np.mean(np.abs(v[idx] - y[idx]))),
                        bias=float(np.mean(v[idx] - y[idx]))) for k, v in p.items()}
    result = dict(n_samples=len(y), n_groups=len(unique), metrics=metrics(np.arange(len(y))),
                  weighting="observation-weighted; paired group bootstrap",
                  evidence_warning="few independent groups" if len(unique) < 10 else None,
                  rmse_improvement_vs_baseline={})
    draws = {k: [] for k in REQUIRED_MODELS if k != "pigment_cell"}
    if len(unique) >= 2:
        rng = np.random.default_rng(seed)
        indices = [np.flatnonzero(g == group) for group in unique]
        for _ in range(bootstrap):
            idx = np.concatenate([indices[i] for i in rng.integers(len(unique), size=len(unique))])
            m = metrics(idx)
            for k in draws:
                draws[k].append(m[k]["rmse"] - m["pigment_cell"]["rmse"])
    for k, values in draws.items():
        result["rmse_improvement_vs_baseline"][k] = dict(
            estimate=result["metrics"][k]["rmse"] - result["metrics"]["pigment_cell"]["rmse"],
            ci95=list(map(float, np.quantile(values, [0.025, 0.975]))) if values else None)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("predictions", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = json.loads(args.manifest.read_text())
    required = {"quantity", "units", "observation_source", "training_ids", "training_groups",
                "model_provenance", "calibration_frozen_before_test"}
    if not required <= manifest.keys() or manifest["calibration_frozen_before_test"] is not True:
        raise ValueError("manifest must document quantity/units/source/provenance and frozen calibration")
    if set(manifest["model_provenance"]) != set(REQUIRED_MODELS):
        raise ValueError("document code/config/calibration provenance for all four models")
    if any(not manifest[k] for k in ("quantity", "units", "observation_source")):
        raise ValueError("quantity, units and source must be nonempty")
    if any(not manifest["model_provenance"][k] for k in REQUIRED_MODELS):
        raise ValueError("model provenance must be nonempty")
    with args.predictions.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    audit_split(manifest["training_ids"], [r["sample_id"] for r in rows],
                manifest["training_groups"], [r["group_id"] for r in rows])
    result = score([float(r["observed"]) for r in rows],
                   {k: [float(r[k]) for r in rows] for k in REQUIRED_MODELS},
                   [r["group_id"] for r in rows])
    result.update(manifest=manifest, predictions_sha256=hashlib.sha256(args.predictions.read_bytes()).hexdigest(),
                  manifest_sha256=hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
                  claim_status="held-out score only; independence and actual-melt validity require source audit")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
