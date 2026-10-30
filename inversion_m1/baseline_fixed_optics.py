#!/usr/bin/env python
"""
Baseline A: fixed-optics conventional inversion of spectral albedo (milestone 1).

Per plot, bounded least squares on the 400-1300 nm hemispherical albedo with BioSNICAR's shipped empirical
optics for glacier algae, Greenland dust and black carbon (a proxy for all dark particles) and a bubbly-ice
column of continuous optical radius. Nothing is fitted across plots; the only hyper-parameter (error scale) is
set from TRAINING days of a leave-one-day-out split.

Stages (each resumable; per-plot records are appended to a JSONL file keyed by the configuration hash):
  fit        per-plot fits: full-window optimum, profile of the algal term, 400-1000 nm fit with the
             1000-1300 nm window held out
  evaluate   leave-one-day-out error scale, intervals, component recovery, broadband albedo, ambiguity metrics
  transfer   the same model on southern-Greenland 2021 HCRF (extra free anisotropy factor k); a LABELLED check on
             reused data, not an untouched test

    python common/run_budgeted.py --background --threads 1 --mem-mb 1500 --name m1_fit -- \
        python inversion_m1/baseline_fixed_optics.py fit
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "phase2"))

import operators as O  # noqa: E402
import splits as SP  # noqa: E402

SCR = "/tmp/claude-0/-home-user-algae/32b7214c-bd8e-5e37-94c0-4ae97c3c08dd/scratchpad"
BIOSNICAR = os.environ.get("BIOSNICAR_PATH", os.path.join(SCR, "ext", "biosnicar-py"))
CONFIG = os.path.join(HERE, "config_baseline_A.json")
OUT_DEFAULT = os.path.join(REPO, "records", "inversion_m1", "baseline_A")


def sha256_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def load_config(path=CONFIG):
    cfg = json.load(open(path))
    cfg["_sha256"] = sha256_file(path)
    return cfg


# ----------------------------------------------------------------------------------------------- data
def s6_plots(cfg):
    """The 47 co-located development plots: albedo spectrum, count, sampling day, noon solar zenith."""
    import empirical_data as ED
    alb = pd.read_csv(os.path.join(BIOSNICAR, "data", "additional_data", "Albedo_master.csv"))
    ct = pd.read_csv(os.path.join(REPO, "data", "empirical", "cook2020_archive_cell_counts.csv"))
    groups = SP.load_groups(os.path.join(REPO, "records", "inversion_m1", "sample_groups.csv"))
    SP.validate_roles(groups)
    dev = groups[(groups.dataset == "S6_2017_albedo") & (groups.role == "development")]
    wl = alb["Wavelength"].to_numpy(float)
    rows = []
    for r in ct.itertuples():
        if r.sample not in set(dev.sample_id):
            continue
        day, month = (int(x) for x in r.sample.split("_")[:2])
        doy = pd.Timestamp(year=2017, month=month, day=day).dayofyear
        rows.append(dict(sample=r.sample, group=f"S6:2017-{month:02d}-{day:02d}", cells=float(r.cells_per_ml),
                         cells_counted=int(r.cells_counted), sza=ED.solar_zenith_noon(doy), spec=alb[r.sample].to_numpy(float)))
    return wl, rows


def bands_for(F, cfg, window):
    lo, hi = window
    return O.window_mask(F.wl_nm, lo, hi)


def observe(wl, spec, F, mask):
    obs = O.band_average(wl, spec, F.wl_nm[mask], 10.0)
    return obs


# ----------------------------------------------------------------------------------------------- fitting
def make_problem(F, params, sza, mask, obs, fixed=None, priors=None, with_k=False, k_prior=None):
    """Residual function in unit-cube coordinates. `fixed` maps parameter index -> fixed physical value."""
    names = [p.name for p in params]
    free = [i for i in range(len(params)) if not (fixed and i in fixed)]

    def theta_of(u):
        th = np.empty(len(params))
        for ui, i in zip(u, free):
            th[i] = params[i].from_unit(ui)
        for i, v in (fixed or {}).items():
            th[i] = v
        return th

    def model(th):
        a = F.albedo(th[0], th[1], th[2], th[3], sza)[mask]
        return th[4] * a if with_k else a

    def resid(u):
        th = theta_of(u)
        r = model(th) - obs
        extra = []
        for i, (mu, sd) in (priors or {}).items():
            extra.append((th[i] - mu) / sd)
        if with_k and k_prior is not None:
            extra.append((th[4] - k_prior[0]) / k_prior[1])
        return np.concatenate([r, np.asarray(extra)]) if extra else r

    return theta_of, model, resid, free


def fit_problem(resid, free, params, cfg, seed, u_start=None):
    """Sobol pre-search then local bounded least squares from the best distinct starts; returns (res, u)."""
    from scipy.stats import qmc
    d = len(free)
    if u_start is not None:
        starts = [np.asarray(u_start, dtype=float)]
    else:
        m = int(np.log2(cfg["fit"]["n_sobol_candidates"]))
        cand = qmc.Sobol(d=d, scramble=True, seed=seed).random_base2(m=m)
        sse = np.array([np.sum(resid(c) ** 2) for c in cand])
        order = np.argsort(sse)
        starts = []
        for i in order:
            if all(np.max(np.abs(cand[i] - s)) > 0.15 for s in starts):
                starts.append(cand[i])
            if len(starts) >= cfg["fit"]["n_local_starts"]:
                break
    best = None
    for s in starts:
        r = least_squares(resid, s, bounds=(0.0, 1.0), method="trf", diff_step=1e-3, max_nfev=cfg["fit"]["max_nfev"],
                          xtol=1e-6, ftol=1e-8)
        if best is None or r.cost < best.cost:
            best = r
    return best


def jtj_theta(res, free, params):
    w = np.array([params[i].hi - params[i].lo for i in free])
    J = res.jac[:, :len(free)] / w[None, :]
    return (J.T @ J), w


def fit_plot(F, cfg, params, wl, rec, seed):
    t0 = time.time()
    fw = cfg["observation"]["fit_window_nm"]
    hw = cfg["observation"]["holdout_check"]
    m_full = bands_for(F, cfg, fw)
    obs_full = observe(wl, rec["spec"], F, m_full)
    ok = np.isfinite(obs_full)
    if ok.mean() < 0.9:
        return dict(sample=rec["sample"], status="rejected", reason=f"only {ok.mean():.0%} of bands finite")
    O.check_albedo_bounds(obs_full[ok], name=rec["sample"])
    m_use = m_full.copy()
    m_use[np.flatnonzero(m_full)[~ok]] = False
    obs = obs_full[ok]
    theta_of, model, resid, free = make_problem(F, params, rec["sza"], m_use, obs)
    res = fit_problem(resid, free, params, cfg, seed)
    theta = theta_of(res.x)
    jt, widths = jtj_theta(res, free, params)
    out = dict(sample=rec["sample"], group=rec["group"], status="ok", sza=rec["sza"], n_bands=int(m_use.sum()),
               theta=theta.tolist(), sse=float(2 * res.cost), nfev=int(res.nfev), at_bound=[bool(params[i].at_bound(theta[i])) for i in range(len(params))],
               jtj=jt.tolist(), cond=float(np.linalg.cond(jt)), resid_full=(model(theta) - obs).tolist())
    # profile of the algal term
    prof = []
    u_best = res.x.copy()
    for g in cfg["fit"]["profile"]["grid"]:
        gl = params[0].lo <= g <= params[0].hi
        if not gl:
            prof.append(None)
            continue
        th_f, mod_f, res_f, free_f = make_problem(F, params, rec["sza"], m_use, obs, fixed={0: g})
        u0 = np.array([params[i].to_unit(theta[i]) for i in free_f])
        r_f = fit_problem(res_f, free_f, params, cfg, seed + 1, u_start=u0)
        prof.append(float(2 * r_f.cost))
    out["profile_grid"] = list(cfg["fit"]["profile"]["grid"])
    out["profile_sse"] = prof
    # hold-out window: fit 400-1000, predict 1000-1300
    m_a = bands_for(F, cfg, hw["fit_window_nm"])
    m_a_use = m_a & m_use
    obs_a = observe(wl, rec["spec"], F, m_a_use)
    th_a, mod_a, res_a, free_a = make_problem(F, params, rec["sza"], m_a_use, obs_a)
    r_a = fit_problem(res_a, free_a, params, cfg, seed + 2)
    theta_a = th_a(r_a.x)
    p_lo, p_hi = hw["predict_window_nm"]
    m_p = O.window_mask(F.wl_nm, p_lo + 1e-9, p_hi) & m_use
    obs_p = observe(wl, rec["spec"], F, m_p)
    pred = F.albedo(*theta_a[:4], rec["sza"])[m_p]
    out.update(theta_partial=theta_a.tolist(), sse_partial=float(2 * r_a.cost), holdout_resid=(pred - obs_p).tolist(),
               holdout_bands_nm=F.wl_nm[m_p].tolist())
    # broadband albedo: observed (H1 convention) vs predicted by the full-window fit
    spec_full = F.albedo(*theta[:4], rec["sza"])
    fl = F.flux(rec["sza"])
    out["bba_obs"] = O.broadband(wl, rec["spec"], F.wl_nm, fl)
    out["bba_pred"] = O.broadband(F.wl_nm, spec_full, F.wl_nm, fl)
    out["seconds"] = round(time.time() - t0, 1)
    return out


def cmd_fit(a):
    from forward import FixedOpticsForward
    cfg = load_config()
    params = O.make_params(cfg["parameters"])
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, "per_plot_fits.jsonl")
    done = {}
    if os.path.isfile(path):
        for ln in open(path):
            d = json.loads(ln)
            if d.get("config_sha256") == cfg["_sha256"]:
                done[d["sample"]] = d
    F = FixedOpticsForward(BIOSNICAR, cfg)
    wl, rows = s6_plots(cfg)
    todo = [r for r in rows if r["sample"] not in done]
    print(f"{len(rows)} plots; {len(done)} done; {len(todo)} to fit", flush=True)
    with open(path, "a") as fh:
        for i, rec in enumerate(todo):
            if a.limit is not None and i >= a.limit:
                break
            out = fit_plot(F, cfg, params, wl, rec, cfg["fit"]["seed"] + i * 7)
            out["config_sha256"] = cfg["_sha256"]
            fh.write(json.dumps(out) + "\n")
            fh.flush()
            print(f"  {rec['sample']:10s} {out.get('status')} cells10^{out['theta'][0]:.2f} (count {rec['cells']:.0f}) "
                  f"rmse {np.sqrt(np.mean(np.square(out['resid_full']))):.4f} {out['seconds']}s", flush=True)


# --------------------------------------------------------------------------------------------- evaluation
def read_fits(out, cfg):
    recs = {}
    for ln in open(os.path.join(out, "per_plot_fits.jsonl")):
        d = json.loads(ln)
        if d.get("config_sha256") == cfg["_sha256"] and d.get("status") == "ok":
            recs[d["sample"]] = d
    return recs


def error_scale_lodo(fits, order, groups, sigma_obs):
    """Leave-one-day-out error scale and residual correlation. For held-out day d, ONLY other days' residuals
    are read (TrainOnly enforces it): sigma_total = max(sigma_obs, RMS of their held-out-window residuals);
    rho = pooled lag-1 autocorrelation of their full-window residuals."""
    hold = [np.asarray(fits[s]["holdout_resid"]) for s in order]
    full = [np.asarray(fits[s]["resid_full"]) for s in order]
    out = {}
    for day, tr, te in SP.leave_one_group_out(groups):
        view_h, view_f = SP.TrainOnly(hold, tr), SP.TrainOnly(full, tr)
        rms = float(np.sqrt(np.mean(np.concatenate([view_h[int(i)] for i in tr]) ** 2)))
        num = den = 0.0
        for i in tr:
            r = view_f[int(i)] - np.mean(view_f[int(i)])
            num += float(np.sum(r[:-1] * r[1:]))
            den += float(np.sum(r * r))
        rho = float(np.clip(num / den, 0.0, 0.95))
        kappa = float(np.clip(np.sqrt((1 + rho) / (1 - rho)), 1.0, 6.0))
        out[day] = dict(sigma_total=max(sigma_obs, rms), rho=rho, kappa=kappa, n_train_plots=int(len(tr)))
    return out


def profile_interval(grid, sse, sse_min, scale):
    """95 % (1-d) profile interval of log10 cells from SSE values; open ends are reported as such."""
    g = np.array([x for x, s in zip(grid, sse) if s is not None], dtype=float)
    d = (np.array([s for s in sse if s is not None]) - sse_min) / scale
    inside = d <= 3.84
    if not inside.any():
        return None, None, True, True
    lo_i, hi_i = np.flatnonzero(inside)[[0, -1]]
    lo = g[lo_i] if lo_i == 0 else float(np.interp(3.84, [d[lo_i], d[lo_i - 1]], [g[lo_i], g[lo_i - 1]]))
    hi = g[hi_i] if hi_i == len(g) - 1 else float(np.interp(3.84, [d[hi_i], d[hi_i + 1]], [g[hi_i], g[hi_i + 1]]))
    return float(lo), float(hi), bool(lo_i == 0), bool(hi_i == len(g) - 1)


def cmd_evaluate(a):
    cfg = load_config()
    params = O.make_params(cfg["parameters"])
    fits = read_fits(a.out, cfg)
    wl, rows = s6_plots(cfg)
    order = [r["sample"] for r in rows if r["sample"] in fits]
    meta = {r["sample"]: r for r in rows}
    groups = [meta[s]["group"] for s in order]
    sigma_obs = cfg["observation"]["error_model"]["sigma_obs"]
    sc = error_scale_lodo(fits, order, groups, sigma_obs)
    recs = []
    for s in order:
        f, m = fits[s], meta[s]
        day = m["group"]
        sig, kap = sc[day]["sigma_total"], sc[day]["kappa"]
        scale = (sig * kap) ** 2
        lo, hi, lo_open, hi_open = profile_interval(f["profile_grid"], f["profile_sse"], f["sse"], scale)
        th = np.array(f["theta"])
        cov = np.linalg.pinv(np.array(f["jtj"])) * scale
        sd = np.sqrt(np.clip(np.diag(cov), 0, None))
        corr = cov / np.outer(np.maximum(sd, 1e-12), np.maximum(sd, 1e-12))
        n_zero = m["cells"] == 0
        n_cnt = max(m["cells_counted"], 1)
        meas = np.log10(m["cells"]) if m["cells"] > 0 else np.nan
        poisson_dex = (1.0 / np.sqrt(n_cnt)) / np.log(10.0)
        hold = np.asarray(f["holdout_resid"])
        recs.append(dict(
            sample=s, group=day, cells_measured=m["cells"], cells_counted=m["cells_counted"], log10_measured=meas,
            log10_retrieved=th[0], log10_dust=th[1], log10_dark=th[2], ln_radius=th[3],
            profile_lo=lo, profile_hi=hi, profile_lo_open=lo_open, profile_hi_open=hi_open,
            laplace_sd_log10_cells=sd[0], laplace_sd_dust=sd[1], laplace_sd_dark=sd[2],
            corr_cells_dust=corr[0, 1], corr_cells_dark=corr[0, 2], corr_cells_radius=corr[0, 3], corr_dust_dark=corr[1, 2],
            dust_at_bound=f["at_bound"][1], dark_at_bound=f["at_bound"][2], cells_at_bound=f["at_bound"][0],
            radius_at_bound=f["at_bound"][3], cond=f["cond"],
            rmse_full=float(np.sqrt(np.mean(np.square(f["resid_full"])))), mae_full=float(np.mean(np.abs(f["resid_full"]))),
            mae_holdout=float(np.mean(np.abs(hold))), bias_holdout=float(np.mean(hold)),
            holdout_z_coverage95=float(np.mean(np.abs(hold) <= 1.96 * sig * kap * 0 + 1.96 * sig)),
            sigma_total=sig, kappa=kap, poisson_sd_dex=poisson_dex,
            bba_obs=f["bba_obs"], bba_pred=f["bba_pred"], bba_err=f["bba_pred"] - f["bba_obs"]))
    D = pd.DataFrame(recs)
    pos = D[D.cells_measured > 0].copy()
    cover = ((pos.log10_measured >= pos.profile_lo - 2 * pos.poisson_sd_dex) & (pos.log10_measured <= pos.profile_hi + 2 * pos.poisson_sd_dex))
    err = pos.log10_retrieved - pos.log10_measured
    sp = stats.spearmanr(pos.log10_retrieved, pos.log10_measured)
    summary = dict(
        config_sha256=cfg["_sha256"], n_plots=int(len(D)), n_days=int(D.group.nunique()), n_positive_counts=int(len(pos)),
        spectral_fit_full_window=dict(median_rmse=float(D.rmse_full.median()), iqr=[float(D.rmse_full.quantile(.25)), float(D.rmse_full.quantile(.75))],
                                      median_mae=float(D.mae_full.median())),
        holdout_window_1000_1300=dict(median_mae=float(D.mae_holdout.median()), mean_bias=float(D.bias_holdout.mean()),
                                      coverage_of_95pct_band_leave_one_day_out=float(D.holdout_z_coverage95.mean())),
        algal_recovery_vs_counts=dict(spearman=float(sp.statistic), spearman_p=float(sp.pvalue), bias_dex=float(err.mean()),
                                      rmse_dex=float(np.sqrt(np.mean(err ** 2))), median_abs_err_dex=float(np.median(np.abs(err))),
                                      interval_coverage_95=float(cover.mean()), profile_interval_median_width_dex=float((pos.profile_hi - pos.profile_lo).median()),
                                      fraction_interval_open_at_either_end=float((pos.profile_lo_open | pos.profile_hi_open).mean()),
                                      fraction_interval_wider_than_1_dex=float(((pos.profile_hi - pos.profile_lo) > 1.0).mean())),
        zero_count_plots=dict(n=int((D.cells_measured == 0).sum()), retrieved_log10_median=float(D[D.cells_measured == 0].log10_retrieved.median()) if (D.cells_measured == 0).any() else None),
        broadband_albedo=dict(bias=float(D.bba_err.mean()), mae=float(D.bba_err.abs().mean()), n=int(len(D))),
        ambiguity=dict(fraction_dust_at_bound=float(D.dust_at_bound.mean()), fraction_dark_particles_at_bound=float(D.dark_at_bound.mean()),
                       fraction_cells_at_bound=float(D.cells_at_bound.mean()), fraction_radius_at_bound=float(D.radius_at_bound.mean()),
                       median_abs_corr_cells_dust=float(D.corr_cells_dust.abs().median()), median_abs_corr_cells_dark=float(D.corr_cells_dark.abs().median()),
                       median_abs_corr_dust_dark=float(D.corr_dust_dark.abs().median()), median_condition_number=float(D.cond.median())),
        error_scale_lodo=dict(sigma_total_range=[min(v["sigma_total"] for v in sc.values()), max(v["sigma_total"] for v in sc.values())],
                              kappa_range=[min(v["kappa"] for v in sc.values()), max(v["kappa"] for v in sc.values())]))
    # published reference: Cook et al. 2020 retrievals shipped with the pinned checkout (29 plots)
    meta_p = pd.read_csv(os.path.join(REPO, "data", "empirical", "cook2020_field_metadata.csv")).set_index("filename")
    both = D.set_index("sample").join(meta_p[["algae_cells_inv_model_S2", "measured_cells"]], how="inner")
    both = both[(both.algae_cells_inv_model_S2 > 0) & (both.cells_measured > 0) & (both.cells_measured.notna())]
    if len(both) > 5:
        pub = np.log10(both.algae_cells_inv_model_S2)
        mine = both.log10_retrieved
        tru = both.log10_measured
        summary["vs_published_cook2020_retrievals"] = dict(
            n=int(len(both)), spearman_mine_vs_published=float(stats.spearmanr(mine, pub).statistic),
            spearman_published_vs_counts=float(stats.spearmanr(pub, tru).statistic), spearman_mine_vs_counts=float(stats.spearmanr(mine, tru).statistic),
            rmse_dex_published_vs_counts=float(np.sqrt(np.mean((pub - tru) ** 2))), rmse_dex_mine_vs_counts=float(np.sqrt(np.mean((mine - tru) ** 2))),
            note="published retrievals are model output of an unverified procedure: a reference, never a truth target")
    summary["provenance"] = dict(generated_utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                                 numpy=np.__version__, scipy=__import__("scipy").__version__, biosnicar_commit=_git(BIOSNICAR),
                                 code_sha256=sha256_file(__file__), operators_sha256=sha256_file(os.path.join(HERE, "operators.py")),
                                 sample_groups_roles_hash=SP.roles_hash(SP.load_groups(os.path.join(REPO, "records", "inversion_m1", "sample_groups.csv"))))
    D.to_csv(os.path.join(a.out, "per_plot_results.csv"), index=False)
    json.dump(dict(per_day_error_scale=sc), open(os.path.join(a.out, "error_scale_lodo.json"), "w"), indent=1)
    json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


def _git(p):
    import subprocess
    try:
        return subprocess.run(["git", "-C", p, "rev-parse", "HEAD"], capture_output=True, text=True, timeout=20).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit")
    f.add_argument("--out", default=OUT_DEFAULT)
    f.add_argument("--limit", type=int, default=None)
    e = sub.add_parser("evaluate")
    e.add_argument("--out", default=OUT_DEFAULT)
    a = p.parse_args(argv)
    {"fit": cmd_fit, "evaluate": cmd_evaluate}[a.cmd](a)


if __name__ == "__main__":
    main()
