#!/usr/bin/env python
"""
Controlled forward-optics diagnostics (Phase 5 item 3): WHY do the physics models over-predict the
measured broadband HCRF of the S6 2017 plots by +0.10 to +0.26 (records/heldout_v2)?

Optics-only with KNOWN state: abundance = the measured count of each plot; every nuisance (bubble radius,
crust density, dust, species mix, illumination) is held at a reference value fixed A PRIORI from the
literature priors (phase4/priors.py), never fitted to the scored spectra. Factors are then changed ONE AT
A TIME. Full retrieval-plus-forward prediction is a different experiment (records/heldout_v2, H1).

Groups (the matrix in docs/forward_diagnostics.md):
  A geometry      measured HCRF vs measured hemispherical albedo of the SAME plot (observation only):
                  k_obs = HCRF/albedo, broadband and per S2 band, vs the retrieval prior k ~ N(0.90, 0.175);
                  exact split of the HCRF error  k_p A_m - k_o A_o = (k_p - k_o) A_m + k_o (A_m - A_o)
                  into a geometry term and an albedo term.
  B ice/impurity  bubble radius (prior 10/90 %), crust density (Cooper et al. 2018 range), crust depth,
                  dust (none / prior 90 %), species mix (f_n 0 / 1), algae in a 1 mm film, granular ice.
  C biology       no algae, Tier A, measured in vivo MACs, TD-DFT tiers C/D, iid calibration.
  D illumination  SZA +/- 5 deg (measurements at noon +/- 2 h), fully diffuse, mid-latitude-summer
                  spectrum, flat (unweighted) broadband.
Not representable in the current bridge (reported as NOT TESTED, not silently skipped): liquid water films
and ponds, surface roughness, snow patches, vertical dust profiles, cloud spectra other than diffuse.

Density diagnostics: prior-ensemble predictive (common random numbers over the nuisance prior) of the
broadband albedo for each optics model; bias vs spread (z = err / sd), PIT and 90 % coverage separate a
biased predictive from a too-narrow one.

    python common/run_budgeted.py --name forward_diag --mem-mb 1500 --threads 1 -- \
        python3 phase4/forward_diagnostics.py --outdir phase4/results/forward_diagnostics
"""

from __future__ import annotations

import argparse
import copy
import dataclasses
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import empirical_data as ED  # noqa: E402

S2 = ("B2", "B3", "B4", "B8")
OPTICS = ("no_algae", "tierA_empirical", "measured_mac_C", "tddft_C", "tddft_D", "tddft_D_iid")
FOCUS = ("tierA_empirical", "tddft_D", "measured_mac_C")          # one-factor matrix
N_ENSEMBLE = 40


# --------------------------------------------------------------------------- observations
def observations(biosnicar_root, runner):
    """Plots with a count, an HCRF and an albedo spectrum; broadband (300-2500 nm, clear-sky irradiance
    at the plot SZA, 300-350 nm held at the 350 nm value) and S2-band values of both quantities."""
    tab, spectra = ED.field_samples(include_stibal=False)
    tab = tab[tab.dataset == "s6_2017"].reset_index(drop=True)
    alb = pd.read_csv(os.path.join(biosnicar_root, "data", "additional_data", "Albedo_master.csv"))
    hc = spectra["s6_2017"]
    srf = pd.read_csv(ED.path("esa_s2_srf_TN-15-0007_v4.0.csv")).set_index("wavelength_nm")
    grid = np.arange(300.0, 2500.0 + 0.5, 1.0)
    rows, excluded, spec_rows = [], [], {}
    for r in tab.itertuples():
        if r.sample not in alb or r.sample not in hc:
            excluded.append(dict(sample=r.sample, reason="no albedo spectrum" if r.sample not in alb else "no HCRF"))
            continue
        a = np.interp(grid, alb.Wavelength, alb[r.sample].to_numpy(float), left=np.nan, right=np.nan)
        h = np.interp(grid, hc.wavelength_nm, hc[r.sample].to_numpy(float), left=np.nan, right=np.nan)
        meas = grid >= 350
        bad = [n for n, v in (("albedo", a), ("HCRF", h)) if np.isfinite(v[meas]).mean() < 0.95]
        if bad:
            excluded.append(dict(sample=r.sample, reason=f">5% missing values in 350-2500 nm ({', '.join(bad)})"))
            continue
        for v in (a, h):                                    # hold the first valid value below 350 nm; fill gaps
            ok = np.isfinite(v)
            v[~ok] = np.interp(grid[~ok], grid[ok], v[ok])
        vis = (grid >= 400) & (grid <= 1300)
        if a[vis].max() > 1.05 or a[vis].min() < 0:
            excluded.append(dict(sample=r.sample, reason="albedo outside [0, 1.05] in 400-1300 nm"))
            continue
        flx = np.interp(grid, runner.wvl_um * 1000.0, runner.illumination(r.sza).flx_slr)
        row = dict(sample=r.sample, day="_".join(r.sample.split("_")[:2]), cells=float(r.cells),
                   sza=float(r.sza), alb_bb=float(np.sum(a * flx) / flx.sum()),
                   hcrf_bb=float(np.sum(h * flx) / flx.sum()), alb_flat=float(a.mean()))
        for b in S2:
            s = srf[f"S2A_{b}"].reindex(grid).fillna(0).to_numpy() * flx
            row[f"alb_{b}"] = float(np.sum(a * s) / s.sum())
            row[f"hcrf_{b}"] = float(np.sum(h * s) / s.sum())
            row[f"k_{b}"] = row[f"hcrf_{b}"] / row[f"alb_{b}"]
        row["k_bb"] = row["hcrf_bb"] / row["alb_bb"]
        rows.append(row)
        spec_rows[r.sample] = a
    return pd.DataFrame(rows), pd.DataFrame(excluded), grid, spec_rows


# --------------------------------------------------------------------------- forward model
@dataclasses.dataclass
class State:
    r_um: float
    rho_top: float = 450.0
    rho_bottom: float = 690.0
    dz_top: float = 0.02
    dust_ppb: float = 0.0
    f_n: float = 0.5
    film_dz: float | None = None
    mode: str = "bubbly"
    sza_offset: float = 0.0


CAL_CACHE = os.path.join(HERE, "results", "cache", "calibrations")


def _with_disk_calibration_cache(fn):
    """Persist tddft_calibration's in-process cache across runs. Entries are keyed by that module's own
    content fingerprint (spectrum, data, code), so a stale entry can never be reused; physics code is
    not modified."""
    import pickle
    import tddft_calibration as TC
    os.makedirs(CAL_CACHE, exist_ok=True)
    for f in os.listdir(CAL_CACHE):
        k = f[:-4]
        if f.endswith(".pkl") and k not in TC._CACHE:
            try:
                TC._CACHE[k] = pickle.load(open(os.path.join(CAL_CACHE, f), "rb"))
            except Exception:  # noqa: BLE001 - unreadable cache entry: recompute
                pass
    out = fn()
    for k, v in TC._CACHE.items():
        dst = os.path.join(CAL_CACHE, f"{k}.pkl")
        if not os.path.exists(dst):
            tmp = dst + f".tmp{os.getpid()}"
            pickle.dump(v, open(tmp, "wb"))
            os.replace(tmp, dst)
    return out


class Forward:
    """Direct BioSNICAR (no emulator) for one optics model; illumination variants via separate runners."""

    def __init__(self, optics, biosnicar, phase1_l2):
        import emulator as E
        import heldout as HO
        self.optics = optics
        key = "tierA_empirical" if optics == "no_algae" else optics
        cfg = E.EmulatorConfig(sza=47, spacecraft="S2A", photosynthetic=True, **HO.OPTICS[key])
        self.b = _with_disk_calibration_cache(lambda: E._Builder(cfg, phase1_l2, False, biosnicar))
        self.bb = self.b.bb
        self.runners = {}
        self.root = self.bb.locate_biosnicar(biosnicar)

    def runner(self, incoming=3, direct=1):
        k = (incoming, direct)
        if k not in self.runners:
            self.runners[k] = self.bb.BioSNICARRunner(self.root, incoming=incoming, direct=direct)
        return self.runners[k]

    def albedo(self, cells, sza, st: State, incoming=3, direct=1):
        """Spectral albedo (480,), normalised irradiance (480,) at the plot."""
        run = self.runner(incoming, direct)
        spec = self.bb.IceSpec(st.r_um, st.rho_top, st.rho_bottom, dz_top=st.dz_top, mode=st.mode,
                               film_dz=st.film_dz, film_only=True)
        imps = []
        if self.optics != "no_algae" and cells > 0:
            if self.b.cfg.model == "ours":
                imps = [(copy.copy(self.b.imps["nordenskioeldii"]), cells * st.f_n),
                        (copy.copy(self.b.imps["alaskanum"]), cells * (1 - st.f_n))]
            else:
                imps = [(copy.copy(self.b.imps["tierA"]), cells)]
        if st.dust_ppb > 0:
            imps.append((copy.copy(self.b.dust), st.dust_ppb))
        alb, flx, _ = run.run_multi(spec, sza + st.sza_offset, imps)
        return alb, flx

    def summary(self, alb, flx, flat=False):
        run = self.runner()
        if flat:
            return float(alb[run.band].mean())
        out = dict(bba=run.broadband(alb, flx))
        out.update({f"m_{b}": v for b, v in zip(S2, self.b._bands(alb, flx))})
        return out


def reference(prior):
    return State(r_um=float(np.exp(prior.mu_lnr)), dust_ppb=float(np.exp(prior.mu_lndust)),
                 f_n=float(prior.f_alpha / (prior.f_alpha + prior.f_beta)))


def variants(prior, ref):
    q = stats.norm.ppf
    v = {"reference": (ref, {})}
    v["B: radius prior 10%"] = (dataclasses.replace(ref, r_um=float(np.exp(prior.mu_lnr + q(0.1) * prior.sd_lnr))), {})
    v["B: radius prior 90%"] = (dataclasses.replace(ref, r_um=float(np.exp(prior.mu_lnr + q(0.9) * prior.sd_lnr))), {})
    lo, hi = ED.ICE_DENSITY["weathering_crust"]["lo"], ED.ICE_DENSITY["weathering_crust"]["hi"]
    v[f"B: crust density {lo:.0f}"] = (dataclasses.replace(ref, rho_top=lo), {})
    v[f"B: crust density {hi:.0f}"] = (dataclasses.replace(ref, rho_top=hi), {})
    v["B: crust depth 1 cm"] = (dataclasses.replace(ref, dz_top=0.01), {})
    v["B: crust depth 5 cm"] = (dataclasses.replace(ref, dz_top=0.05), {})
    v["B: no dust"] = (dataclasses.replace(ref, dust_ppb=0.0), {})
    v["B: dust prior 90%"] = (dataclasses.replace(ref, dust_ppb=float(np.exp(prior.mu_lndust + q(0.9) * prior.sd_lndust))), {})
    v["B: f_n = 0 (all A. alaskanum)"] = (dataclasses.replace(ref, f_n=0.0), {})
    v["B: f_n = 1 (all A. nordenskioeldii)"] = (dataclasses.replace(ref, f_n=1.0), {})
    v["B: algae in 1 mm film"] = (dataclasses.replace(ref, film_dz=0.001), {})
    v["B: granular ice (r 1500 um)"] = (dataclasses.replace(ref, mode="grains", r_um=1500.0), {})
    v["D: SZA -5 deg"] = (dataclasses.replace(ref, sza_offset=-5.0), {})
    v["D: SZA +5 deg"] = (dataclasses.replace(ref, sza_offset=5.0), {})
    v["D: fully diffuse"] = (ref, dict(direct=0))
    v["D: mid-latitude summer spectrum"] = (ref, dict(incoming=1))
    return v


def ensemble_states(prior, ref, n, seed=0):
    """Common random numbers over the nuisance prior: radius, dust, f_n (crust held at reference)."""
    rng = np.random.default_rng(seed)
    r = np.exp(prior.mu_lnr + prior.sd_lnr * rng.standard_normal(n))
    d = np.exp(prior.mu_lndust + prior.sd_lndust * rng.standard_normal(n))
    f = rng.beta(prior.f_alpha, prior.f_beta, n)
    return [dataclasses.replace(ref, r_um=float(a), dust_ppb=float(b), f_n=float(c)) for a, b, c in zip(r, d, f)]


# --------------------------------------------------------------------------- run
def run(biosnicar=None, phase1_l2=None, outdir=None, n_ens=N_ENSEMBLE, optics=OPTICS, focus=FOCUS, limit=None):
    import biosnicar_bridge as bb
    from priors import PriorConfig
    root = bb.locate_biosnicar(biosnicar)
    base = bb.BioSNICARRunner(root)
    obs, excluded, grid, aspec = observations(root, base)
    if limit:
        obs = obs.iloc[:limit].reset_index(drop=True)
    prior = PriorConfig.for_density(690.0)
    ref = reference(prior)
    k_p = float(prior.mu_k)
    rows, resid, ens_rows = [], {}, []
    wl = base.wvl_um * 1000.0
    for opt in optics:
        fw = Forward(opt, biosnicar, phase1_l2)
        vs = variants(prior, ref) if opt in focus else {"reference": (ref, {})}
        for vname, (st, ill) in vs.items():
            for r in obs.itertuples():
                alb, flx = fw.albedo(r.cells, r.sza, st, **ill)
                s = fw.summary(alb, flx)
                A_m, A_o, k_o = s["bba"], r.alb_bb, r.k_bb
                rows.append(dict(optics=opt, variant=vname, group=vname.split(":")[0], sample=r.sample, day=r.day,
                                 cells=r.cells, alb_obs=A_o, alb_mod=A_m, err_alb=A_m - A_o,
                                 err_alb_flat=fw.summary(alb, flx, flat=True) - r.alb_flat,
                                 hcrf_obs=r.hcrf_bb, hcrf_pred=k_p * A_m, err_hcrf=k_p * A_m - r.hcrf_bb,
                                 err_hcrf_geometry=(k_p - k_o) * A_m, err_hcrf_albedo=k_o * (A_m - A_o),
                                 **{f"err_{b}": s[f"m_{b}"] - getattr(r, f"alb_{b}") for b in S2}))
                if vname == "reference":
                    resid.setdefault(opt, []).append(np.interp(grid, wl, alb) - aspec[r.sample])
        if n_ens:
            states = ensemble_states(prior, ref, n_ens)
            for r in obs.itertuples():
                draws = np.array([fw.summary(*fw.albedo(r.cells, r.sza, st))["bba"] for st in states])
                m, sd = draws.mean(), draws.std(ddof=1)
                q05, q95 = np.quantile(draws, [0.05, 0.95])
                ens_rows.append(dict(optics=opt, sample=r.sample, day=r.day, alb_obs=r.alb_bb, ens_mean=m, ens_sd=sd,
                                     err=m - r.alb_bb, z=(m - r.alb_bb) / sd if sd > 0 else np.nan,
                                     pit=float(np.mean(draws <= r.alb_bb)), covered90=float(q05 <= r.alb_bb <= q95),
                                     width90=q95 - q05))
        print(f"[forward-diag] {opt} done", flush=True)
    R = pd.DataFrame(rows)
    EN = pd.DataFrame(ens_rows)

    def agg(g):
        return pd.Series(dict(n=len(g), bias_alb=g.err_alb.mean(), mae_alb=g.err_alb.abs().mean(),
                              rmse_alb=np.sqrt((g.err_alb ** 2).mean()), bias_alb_flat=g.err_alb_flat.mean(),
                              bias_hcrf=g.err_hcrf.mean(), hcrf_geometry=g.err_hcrf_geometry.mean(),
                              hcrf_albedo=g.err_hcrf_albedo.mean(),
                              **{f"bias_{b}": g[f"err_{b}"].mean() for b in S2}))
    S = R.groupby(["optics", "variant"], sort=False).apply(agg, include_groups=False).reset_index()
    geo = dict(n=int(len(obs)), k_prior_mean=k_p, k_prior_sd=float(prior.sd_k),
               k_obs_bb_mean=float(obs.k_bb.mean()), k_obs_bb_sd=float(obs.k_bb.std(ddof=1)),
               k_obs_bb_q=[float(x) for x in obs.k_bb.quantile([0.05, 0.5, 0.95])],
               **{f"k_obs_{b}_mean": float(obs[f"k_{b}"].mean()) for b in S2},
               corr_k_vs_logcells=float(np.corrcoef(obs.k_bb, np.log10(obs.cells.clip(lower=10)))[0, 1]),
               alb_obs_bb_mean=float(obs.alb_bb.mean()), hcrf_obs_bb_mean=float(obs.hcrf_bb.mean()),
               note="k_prior is built from ARF_master = HCRF/albedo of S6 2017 plots incl. 29 of these (deviation D3)")
    ES = (EN.groupby("optics").agg(n=("sample", "size"), bias=("err", "mean"), rmse=("err", lambda e: np.sqrt((e ** 2).mean())),
                                   mean_sd=("ens_sd", "mean"), z_mean=("z", "mean"), z_sd=("z", "std"),
                                   coverage90=("covered90", "mean"), width90=("width90", "mean"),
                                   pit_mean=("pit", "mean")).reset_index() if len(EN) else pd.DataFrame())
    RS = pd.DataFrame({"wavelength_nm": grid, **{f"{o}_median": np.median(v, axis=0) for o, v in resid.items()},
                       **{f"{o}_q25": np.quantile(v, 0.25, axis=0) for o, v in resid.items()},
                       **{f"{o}_q75": np.quantile(v, 0.75, axis=0) for o, v in resid.items()}})
    if outdir:
        import provenance as PV
        os.makedirs(outdir, exist_ok=True)
        PV.atomic_to_csv(obs, os.path.join(outdir, "observations.csv"), index=False, float_format="%.5g")
        PV.atomic_to_csv(excluded, os.path.join(outdir, "excluded.csv"), index=False)
        PV.atomic_to_csv(R, os.path.join(outdir, "per_plot.csv"), index=False, float_format="%.5g")
        PV.atomic_to_csv(S, os.path.join(outdir, "summary.csv"), index=False, float_format="%.5g")
        PV.atomic_to_csv(EN, os.path.join(outdir, "ensemble_per_plot.csv"), index=False, float_format="%.5g")
        PV.atomic_to_csv(ES, os.path.join(outdir, "ensemble_summary.csv"), index=False, float_format="%.5g")
        PV.atomic_to_csv(RS.iloc[::5], os.path.join(outdir, "residual_spectra_5nm.csv"), index=False, float_format="%.5g")
        PV.atomic_write_text(os.path.join(outdir, "settings.json"), json.dumps(dict(
            geometry=geo, reference_state=dataclasses.asdict(ref), n_ensemble=n_ens, optics=list(optics),
            focus=list(focus), prior=dataclasses.asdict(prior), environment=PV.environment(),
            not_tested=["liquid water films/ponds", "surface roughness", "snow patches", "vertical dust profile",
                        "cloud spectra other than fully diffuse", "packaging/scattering variants of the cell model",
                        "Fe increment separately from tier D"]), indent=1, default=float))
    return obs, R, S, EN, ES, RS, geo


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "forward_diagnostics"))
    p.add_argument("--n-ensemble", type=int, default=N_ENSEMBLE)
    p.add_argument("--limit", type=int, default=None, help="first N plots only (smoke test)")
    a = p.parse_args(argv)
    obs, R, S, EN, ES, RS, geo = run(a.biosnicar, a.phase1_l2, a.outdir, a.n_ensemble, limit=a.limit)
    pd.set_option("display.width", 250)
    print(json.dumps(geo, indent=1))
    print(S.round(4).to_string(index=False))
    print(ES.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
