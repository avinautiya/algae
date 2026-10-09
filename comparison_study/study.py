"""Small-data, training-day-tuned ML baselines; positive-count supplementary analysis."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

# Set before numerical imports; threadpool_limits also handles preloaded libraries.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import numpy as np
import pandas as pd
import sklearn
import scipy
from sklearn.compose import TransformedTargetRegressor
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

BANDS = ["B2", "B3", "B4", "B8"]
FOLDS = [("primary", "s6_2017", "sgris_2021"), ("secondary", "sgris_2021", "s6_2017")]
CANDIDATES = {
    "mean": [{}],
    "ridge": [{"alpha": a} for a in (.1, 1., 10., 100.)],
    "random_forest": [{"min_samples_leaf": m} for m in (2, 5)],
    "gaussian_process": [{"length_scale": l, "noise": n} for l in (.5, 2.) for n in (.05, .2)],
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(model, params, seed=0):
    if model == "mean":
        est = DummyRegressor(strategy="mean")
    elif model == "ridge":
        est = Ridge(**params)
    elif model == "random_forest":
        est = RandomForestRegressor(n_estimators=100, max_depth=3, random_state=seed, n_jobs=1, **params)
    elif model == "gaussian_process":
        kernel = RBF(params["length_scale"], length_scale_bounds="fixed") + WhiteKernel(
            params["noise"], noise_level_bounds="fixed")
        est = GaussianProcessRegressor(kernel=kernel, optimizer=None, alpha=1e-8, random_state=seed)
    else:
        raise ValueError("unknown model")
    return TransformedTargetRegressor(regressor=make_pipeline(StandardScaler(), est),
                                      transformer=StandardScaler())


def tune(X, y, groups, model, seed=0):
    """Only training data accepted; complete-day inner holdout with fold-local preprocessing."""
    X, y, groups = np.asarray(X, float), np.asarray(y, float), np.asarray(groups, str)
    if y.ndim == 1:
        y = y[:, None]
    if X.ndim != 2 or y.ndim != 2 or len(X) != len(y) or groups.shape != (len(y),):
        raise ValueError("training array shapes disagree")
    if not np.isfinite(X).all() or not np.isfinite(y).all():
        raise ValueError("missing/nonfinite training values")
    days = np.unique(groups)
    if len(days) < 2:
        raise ValueError("at least two training sampling days required")
    scores = []
    for params in CANDIDATES[model]:
        daily = []
        for day in days:
            val = groups == day
            estimator = build(model, params, seed).fit(X[~val], y[~val])
            pred = np.asarray(estimator.predict(X[val])).reshape(y[val].shape)
            daily.append(float(np.mean((pred - y[val]) ** 2)))
        scores.append(dict(params=params, mean_day_mse=float(np.mean(daily)), daily_mse=daily))
    chosen = min(range(len(scores)), key=lambda i: scores[i]["mean_day_mse"])
    estimator = build(model, scores[chosen]["params"], seed).fit(X, y)
    return estimator, dict(chosen=scores[chosen]["params"], candidates=scores, training_days=days.tolist())


def metrics(y, pred):
    y, pred = np.asarray(y, float), np.asarray(pred, float)
    if y.shape != pred.shape or not np.isfinite(y).all() or not np.isfinite(pred).all():
        raise ValueError("finite matched predictions required")
    e = pred - y
    return dict(rmse=float(np.sqrt(np.mean(e ** 2))), mae=float(np.mean(np.abs(e))), bias=float(e.mean()),
                per_output_rmse=np.sqrt(np.mean(e ** 2, axis=0)).tolist(),
                per_output_bias=np.mean(e, axis=0).tolist())


def paired_rmse(y, pred_a, pred_b, groups, n=2000, seed=0):
    """Positive = A has smaller RMSE. Suppress degenerate few-day confidence intervals."""
    y, a, b = (np.asarray(v, float) for v in (y, pred_a, pred_b))
    if y.shape != a.shape or y.shape != b.shape or not all(np.isfinite(v).all() for v in (y, a, b)):
        raise ValueError("paired arrays must match and be finite")
    groups = np.asarray(groups, str)
    if groups.shape != (len(y),):
        raise ValueError("group shape mismatch")
    blocks = np.unique(groups)
    estimate = metrics(y, b)["rmse"] - metrics(y, a)["rmse"]
    out = dict(rmse_improvement=estimate, n_samples=len(y), n_days=len(blocks), ci95=None,
               interpretation="insufficient days" if len(blocks) < 5 else "conditional day-block sampling uncertainty")
    if len(blocks) >= 5:
        rng = np.random.default_rng(seed)
        idx = [np.flatnonzero(groups == day) for day in blocks]
        values = []
        for _ in range(n):
            rows = np.concatenate([idx[i] for i in rng.integers(len(blocks), size=len(blocks))])
            values.append(metrics(y[rows], b[rows])["rmse"] - metrics(y[rows], a[rows])["rmse"])
        out["ci95"] = np.quantile(values, [.025, .975]).tolist()
    return out


def sampling_day(dataset, sample):
    if dataset == "s6_2017":
        day, month, *_ = sample.split("_")
        stamp = pd.Timestamp(year=2017, month=int(month), day=int(day))
    elif dataset == "sgris_2021":
        stamp = pd.to_datetime(sample.split("-")[0], format="%y%m%d")
    else:
        raise ValueError("unsupported dataset")
    return dataset + ":" + stamp.strftime("%Y-%m-%d")


def load_field(root):
    sys.path.insert(0, str(root / "phase4"))
    from field_validation import field_band_reflectance
    df = field_band_reflectance("S2A")
    df = df[df.dataset.isin([f[1] for f in FOLDS])].copy().reset_index(drop=True)
    df["sample_id"] = df.dataset + ":" + df["sample"]
    df["day"] = [sampling_day(d, s) for d, s in zip(df.dataset, df["sample"])]
    df["included_positive_count"] = df.cells > 0
    if df.sample_id.duplicated().any() or not np.isfinite(df[BANDS + ["cells"]].to_numpy()).all():
        raise ValueError("duplicate sample keys or nonfinite source measurements")
    return df


def physics_comparisons(root, fold, samples, observed, groups, ml_pred):
    """Exact sample/target check; positive-count point metrics only, not count-density scores."""
    path = root / "records/heldout_v2/heldout_per_sample.csv"
    if not path.exists():
        return []
    df = pd.read_csv(path)
    out = []
    for method, rows in df[df.fold == fold].groupby("method"):
        rows = rows[np.isfinite(rows.log_b_obs)]
        if rows["sample"].duplicated().any() or set(rows["sample"]) != set(samples):
            raise ValueError(f"physics sample mismatch: {fold}/{method}")
        rows = rows.set_index("sample").loc[list(samples)]
        if not np.allclose(rows.log_b_obs.to_numpy(), observed[:, 0], atol=5e-5, rtol=0):
            raise ValueError("physics observed targets do not match")
        p = rows.median_obs_pred.to_numpy()[:, None]
        for model, values in ml_pred.items():
            out.append(dict(fold=fold, task="inverse", model=model, comparator=method,
                            comparator_prediction="reported physics/baseline median, not mean",
                            **paired_rmse(observed, values, p, groups)))
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    out = args.outdir
    if out.exists() and not args.resume:
        raise FileExistsError("new output directory required unless --resume")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    df = load_field(root)
    data_bytes = df.to_csv(index=False, float_format="%.17g").encode()
    dependency_versions = dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                               scipy=scipy.__version__, sklearn=sklearn.__version__)
    signature = hashlib.sha256(data_bytes + Path(__file__).read_bytes() +
                               json.dumps(dependency_versions, sort_keys=True).encode() +
                               (root / "comparison_study/PROTOCOL.md").read_bytes()).hexdigest()
    old = out / "signature.txt"
    if old.exists() and old.read_text().strip() != signature:
        raise ValueError("resume input/code/dependency signature mismatch")
    old.write_text(signature + "\n")
    df.to_csv(out / "data_inventory.csv", index=False, float_format="%.17g")
    positive = df[df.included_positive_count].copy()
    results, predictions, comparisons, selection = [], [], [], {}
    with threadpool_limits(limits=1):
        for fold, train_site, test_site in FOLDS:
            tr, te = positive[positive.dataset == train_site], positive[positive.dataset == test_site]
            if set(tr.sample_id) & set(te.sample_id) or set(tr.day) & set(te.day):
                raise ValueError("outer sample/day leakage")
            for task in ("inverse", "forward"):
                Xt = tr[BANDS].to_numpy() if task == "inverse" else np.log10(tr.cells.to_numpy())[:, None]
                Xv = te[BANDS].to_numpy() if task == "inverse" else np.log10(te.cells.to_numpy())[:, None]
                yt = np.log10(tr.cells.to_numpy())[:, None] if task == "inverse" else tr[BANDS].to_numpy()
                yv = np.log10(te.cells.to_numpy())[:, None] if task == "inverse" else te[BANDS].to_numpy()
                extrap = np.any((Xv < Xt.min(axis=0)) | (Xv > Xt.max(axis=0)), axis=1)
                pred_by = {}
                for model in CANDIDATES:
                    key = f"{fold}_{task}_{model}"
                    checkpoint = out / (key + ".json")
                    if checkpoint.exists() and args.resume:
                        saved = json.loads(checkpoint.read_text())
                        if saved["signature"] != signature:
                            raise ValueError("checkpoint signature mismatch")
                        pred, sel = np.array(saved["predictions"]), saved["selection"]
                    else:
                        fit, sel = tune(Xt, yt, tr.day.to_numpy(), model)
                        pred = np.asarray(fit.predict(Xv)).reshape(yv.shape)
                        temp = checkpoint.with_suffix(".tmp")
                        temp.write_text(json.dumps(dict(signature=signature, predictions=pred.tolist(), selection=sel),
                                                   indent=2, allow_nan=False))
                        temp.replace(checkpoint)
                    pred_by[model] = pred
                    selection[key] = sel
                    results.append(dict(fold=fold, task=task, model=model, n_train=len(tr), n_test=len(te),
                                        n_train_days=tr.day.nunique(), n_test_days=te.day.nunique(),
                                        n_extrapolation=int(extrap.sum()), units="dex" if task == "inverse" else "HCRF",
                                        **metrics(yv, pred)))
                    for i, row in enumerate(te.itertuples()):
                        for j, output in enumerate(["log10_cells_ml"] if task == "inverse" else BANDS):
                            predictions.append(dict(fold=fold, task=task, model=model, sample_id=row.sample_id,
                                                    day=row.day, output=output, observed=yv[i, j], predicted=pred[i, j],
                                                    extrapolation=bool(extrap[i])))
                    print(key, results[-1]["rmse"], flush=True)
                for model in CANDIDATES:
                    if model != "mean":
                        comparisons.append(dict(fold=fold, task=task, model=model, comparator="training_mean",
                                                **paired_rmse(yv, pred_by[model], pred_by["mean"], te.day.to_numpy())))
                if task == "inverse":
                    comparisons += physics_comparisons(root, fold, te["sample"].to_numpy(), yv, te.day.to_numpy(), pred_by)
    pd.DataFrame(results).to_csv(out / "summary.csv", index=False)
    pd.DataFrame(predictions).to_csv(out / "predictions.csv", index=False)
    (out / "paired_comparisons.json").write_text(json.dumps(comparisons, indent=2, allow_nan=False) + "\n")
    (out / "selection.json").write_text(json.dumps(selection, indent=2, allow_nan=False) + "\n")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    manifest = dict(status="retrospective supplementary positive-count field comparison; not melt/satellite validation",
                    signature=signature, source_commit=commit, code_sha256=sha(__file__),
                    protocol_sha256=sha(root / "comparison_study/PROTOCOL.md"),
                    data_sha256=hashlib.sha256(data_bytes).hexdigest(),
                    empirical_csv_sha256={p.name: sha(p) for p in sorted((root / "data/empirical").glob("*.csv"))},
                    versions=dependency_versions, excluded_zero_sample_ids=df.loc[~df.included_positive_count, "sample_id"].tolist(),
                    physics_comparator_sha256=sha(root / "records/heldout_v2/heldout_per_sample.csv")
                        if (root / "records/heldout_v2/heldout_per_sample.csv").exists() else None,
                    runtime_seconds=time.time() - t0, max_rss_raw=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    max_rss_units="bytes on macOS; KiB on Linux", no_winner_selected=True)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
