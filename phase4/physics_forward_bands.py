#!/usr/bin/env python
"""Matched per-band physics forward predictions on exactly the comparison_study positive-count rows, so physics
four-band HCRF RMSE can be compared with the ML forward RMSE on the same endpoint (Phase 6 item D).

Prediction for a sample with observed log10 abundance:
    E[HCRF_b] = mu_k x E_prior[F_b(log B_obs, f_n, r, dust)]
using the frozen held-out emulators at the sample's SZA node, the FOLD-TRAINED radius prior
(heldout_settings.json, '<fold>/<optics>') and the literature nuisance priors. mu_k is the k-prior mean
(deviation D3: built from S6 2017 ARF, which overlaps the S6 rows of the secondary fold).
This is the physics point prediction (mean); ML models are squared-error point predictors. Both are point
estimates on identical rows and outputs (B2, B3, B4, B8).

    python phase4/physics_forward_bands.py --outdir records/ml_comparison_physics
"""
import argparse
import dataclasses
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "phase2"))
sys.path.insert(0, ROOT)

BANDS = ["B2", "B3", "B4", "B8"]
OPTICS = ("tierA_empirical", "measured_mac_C", "tddft_C", "tddft_D", "tddft_D_iid")
EMU = os.path.join(HERE, "results", "heldout_v2", "cache")


def predict(opt, fold, rows):
    import emulator as E
    from priors import PriorConfig, prior_logpdfs
    s = json.load(open(os.path.join(HERE, "results", "heldout_v2", "heldout_settings.json")))["settings"][f"{fold}/{opt}"]
    pc = dataclasses.replace(PriorConfig.for_density(690.0), mu_lnr=float(np.log(s["radius_median_um"])),
                             sd_lnr=float(s["radius_ln_sd"]))
    out, ems = [], {}
    for r in rows.itertuples():
        z = int(round(r.sza))
        if z not in ems:
            ems[z] = E.Emulator.load(os.path.join(EMU, f"heldout_{opt}_sza{z}.npz")).refine_log_b(0.05)
        em = ems[z]
        lp, sk, mk, _ = prior_logpdfs(em.axes, 1, pc)
        names = em.names[1:]
        w = sum(lp[a][0].reshape([-1 if b == a else 1 for b in names]) for a in names)
        w = np.exp(w - w.max())
        w /= w.sum()
        lb = em.axes["log_b"]
        j = int(np.clip(np.searchsorted(lb, np.log10(r.cells)) - 1, 0, len(lb) - 2))
        t = float(np.clip((np.log10(r.cells) - lb[j]) / (lb[j + 1] - lb[j]), 0, 1))
        F = (1 - t) * em.data["bands"][j] + t * em.data["bands"][j + 1]          # (f_n, r, dust, 4)
        out.append(mk * np.tensordot(w, F, axes=(list(range(w.ndim)), list(range(w.ndim)))))
    return np.array(out)


def main(argv=None):
    from comparison_study.study import paired_rmse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", default=os.path.join(ROOT, "records", "ml_comparison_physics"))
    a = p.parse_args(argv)
    inv = pd.read_csv(os.path.join(ROOT, "records", "ml_comparison", "data_inventory.csv")).set_index("sample_id")
    ml = pd.read_csv(os.path.join(ROOT, "records", "ml_comparison", "predictions.csv"))
    ml = ml[ml.task == "forward"]
    preds, summ, pairs = [], [], []
    for fold in ("primary", "secondary"):
        f = ml[ml.fold == fold]
        ids = list(dict.fromkeys(f.sample_id))
        rows = inv.loc[ids].reset_index()
        obs = rows[BANDS].to_numpy(float)
        groups = np.repeat(rows.day.to_numpy(), 4)
        phys = {o: predict(o, fold, rows) for o in OPTICS}
        mlp = {m: f[f.model == m].pivot(index="sample_id", columns="output", values="predicted").loc[ids, BANDS].to_numpy()
               for m in f.model.unique()}
        check = f[f.model == "mean"].pivot(index="sample_id", columns="output", values="observed").loc[ids, BANDS].to_numpy()
        assert np.allclose(check, obs), "observed targets differ from the ML study's rows"
        allm = {**{f"physics_{o}": v for o, v in phys.items()}, **{f"ml_{m}": v for m, v in mlp.items()}}
        for m, P in allm.items():
            e = P - obs
            summ.append(dict(fold=fold, model=m, n_samples=len(ids), rmse_pooled=float(np.sqrt((e ** 2).mean())),
                             **{f"rmse_{b}": float(np.sqrt((e[:, i] ** 2).mean())) for i, b in enumerate(BANDS)},
                             **{f"bias_{b}": float(e[:, i].mean()) for i, b in enumerate(BANDS)}))
            for k, sid in enumerate(ids):
                for i, b in enumerate(BANDS):
                    preds.append(dict(fold=fold, model=m, sample_id=sid, day=rows.day[k], output=b, observed=obs[k, i], predicted=P[k, i]))
        for o in OPTICS:
            for m in mlp:
                r = paired_rmse(obs.ravel(), phys[o].ravel(), mlp[m].ravel(), groups)
                pairs.append(dict(fold=fold, a=f"physics_{o}", b=f"ml_{m}", **{k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.items()}))
    os.makedirs(a.outdir, exist_ok=True)
    S = pd.DataFrame(summ)
    pd.DataFrame(preds).to_csv(os.path.join(a.outdir, "predictions.csv"), index=False, float_format="%.6g")
    S.to_csv(os.path.join(a.outdir, "summary.csv"), index=False, float_format="%.4g")
    pd.DataFrame(pairs).to_csv(os.path.join(a.outdir, "paired_physics_vs_ml.csv"), index=False, float_format="%.4g")
    json.dump(dict(note="positive-count rows of records/ml_comparison; physics point = k-prior mean x prior-mean bands "
                        "(fold-trained radius prior); D3 applies to the secondary fold", k_prior="N(0.898,0.175) from ARF_master"),
              open(os.path.join(a.outdir, "manifest.json"), "w"), indent=1)
    pd.set_option("display.width", 200)
    print(S[["fold", "model", "n_samples", "rmse_pooled", "rmse_B2", "rmse_B8", "bias_B2", "bias_B8"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
