#!/usr/bin/env python
"""
Paper-figure suite, built ONLY from saved artifacts (no chemistry, no BioSNICAR runs).

    python3 paper_figures/make_figures.py [--only fig2 fig3 ...]

For every figure: <name>.pdf (vector), <name>.png (300 dpi), <name>_data.csv (the plotted values) and an
entry in paper_figures/output/manifest.json (caption, source files with sha256, code commit and script
sha256, configuration, validation status). Figures whose evidence is missing are recorded as BLOCKED.

Model colours (validated categorical palette, light mode; contrast relief = legends/direct labels + CSV):
Tier A (M1) slot 1 blue, TD-DFT D/B3LYP (M3) slot 2 orange, measured MAC (M1b) slot 3 aqua,
CAM-B3LYP slot 4 yellow, simple empirical M2 slot 5 magenta, no algae (M0) neutral grey. Tier C is the
functional's colour with a dashed line (secondary encoding), never a new hue.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import pickle
import subprocess
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(HERE, "output")
sys.path[:0] = [os.path.join(ROOT, "phase2"), os.path.join(ROOT, "phase1"), os.path.join(ROOT, "phase4")]

C = dict(tierA="#2a78d6", tddft_D="#eb6834", measured="#1baf7a", cam="#eda100", m2="#e87ba4", none="#8a8984",
         ink="#1d1d1b", ink2="#52514e", grid="#e4e3df", surface="#fcfcfb")
NAMES = {"tierA_empirical": "Tier A (M1)", "measured_mac_C": "Measured MAC (M1b)", "tddft_D": "TD-DFT D, B3LYP (M3)",
         "tddft_C": "TD-DFT C, B3LYP", "tddft_D_iid": "TD-DFT D, iid calib.", "no_algae": "No algae (M0)"}
COL = {"tierA_empirical": C["tierA"], "measured_mac_C": C["measured"], "tddft_D": C["tddft_D"],
       "tddft_C": C["tddft_D"], "tddft_D_iid": C["tddft_D"], "no_algae": C["none"]}
LS = {"tddft_C": "--", "tddft_D_iid": ":"}

plt.rcParams.update({"font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5, "legend.fontsize": 6.5,
                     "axes.edgecolor": C["ink2"], "axes.labelcolor": C["ink"], "xtick.color": C["ink2"],
                     "ytick.color": C["ink2"], "text.color": C["ink"], "axes.linewidth": 0.6,
                     "axes.grid": True, "grid.color": C["grid"], "grid.linewidth": 0.5, "axes.axisbelow": True,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.5,
                     "legend.frameon": False, "figure.facecolor": "white", "savefig.facecolor": "white",
                     "pdf.fonttype": 42})
MANIFEST, USED = {}, {}


def src(path):
    """Register a source file (sha256) for the current figure and return its absolute path."""
    p = os.path.join(ROOT, path)
    USED[path] = hashlib.sha256(open(p, "rb").read()).hexdigest()
    return p


def commit():
    try:
        return subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except OSError:
        return "unknown"


def save(fig, name, data: pd.DataFrame, caption, status, config=None):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(os.path.join(OUT, f"{name}.pdf"))
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300)
    plt.close(fig)
    data.to_csv(os.path.join(OUT, f"{name}_data.csv"), index=False, float_format="%.6g")
    MANIFEST[name] = dict(caption=caption, validation_status=status, sources=dict(USED), config=config or {},
                          files=[f"{name}.pdf", f"{name}.png", f"{name}_data.csv"], build=build_info())
    USED.clear()


def build_info():
    dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "--", "paper_figures/make_figures.py"],
                           capture_output=True, text=True).stdout.strip() != ""
    return dict(code_commit=commit(), script_sha256=hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                script_uncommitted_changes=dirty, matplotlib=matplotlib.__version__)


def blocked(name, reason):
    MANIFEST[name] = dict(status="BLOCKED", reason=reason, build=build_info())
    USED.clear()


def panel(ax, letter):
    t = ax.get_title()
    if t:
        ax.set_title(rf"$\bf{{{letter}}}$  " + t, fontsize=ax.title.get_fontsize())
    else:
        ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


# =========================================================================== Figure 1
def fig1():
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.set_xlim(0, 100); ax.set_ylim(0, 62); ax.axis("off"); ax.grid(False)
    style = {"computed": ("#fde8dc", C["tddft_D"]), "measured": ("#d9f2e8", C["measured"]),
             "model": ("#e1ecfa", C["tierA"]), "obs": ("#f3f2ee", C["ink2"])}

    def box(x, y, w, h, text, kind, lw=0.8):
        fc, ec = style[kind]
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2", fc=fc, ec=ec, lw=lw))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=5.6)

    def arrow(x0, y0, x1, y1, label=None):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=7, lw=0.7, color=C["ink2"]))
        if label:
            ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 0.8, label, fontsize=5.2, color=C["ink2"], ha="center")
    box(1, 50, 16, 9, "TD-DFT states\n(B3LYP, CAM-B3LYP;\nneutral glucoside, PCM)", "computed")
    box(21, 50, 16, 9, "Broadened MAC\n(Gaussian in energy;\nNapierian, per glucoside)", "computed")
    box(21, 37, 16, 9, "Isolated chromophore\nHPLC shape\n(calibration target)", "measured")
    box(21, 23, 16, 10, "Extract MAC, phenol eq.\n(calibration target)\n+ Fe increment\n(Procházková 2025)", "measured")
    box(44, 37, 14, 11, "Calibration\n(cut posterior)\ndE, w | f, φ\n4000 joint draws", "model")
    box(65, 50, 14, 9, "Cell optics\n(packaging,\nchl/carot., Mie g)", "model")
    box(65, 37, 14, 9, "BioSNICAR column\n(bubbly ice, dust,\nclear-sky SZA 50°)", "model")
    box(84, 37, 15, 9, "Glacier spectral\nalbedo; Δα;\nabsorbed SW", "model")
    box(65, 22, 14, 10, "Observation operator\nalbedo → HCRF = k·α\n(scalar k, D3)\nS2 SRF bands", "model")
    box(84, 22, 15, 10, "Host interface\n(reference SEB,\nModes A/B)", "model")
    box(1, 2, 31, 13, "Validation / development data\nS6 2017 albedo (H1): development-reused\nKAN_L/M radiation + S2 (H2/H4)\nPROMBIO (H5): 2 test station-days", "obs")
    box(35, 2, 29, 13, "Fitted on (never validation)\nWilliamson 2020 MAC, HPLC;\nS6 2017 HCRF (σ, τ, radius prior);\nARF → k prior (overlaps S6, D3)", "obs")
    box(67, 2, 32, 13, "Outcomes (this study)\nH1/H2/H4 not supported; H5 insufficient;\nH3 untested; molecular substitution\nunresolved after calibration", "obs")
    arrow(17.4, 54.5, 20.6, 54.5)
    arrow(37.4, 54.5, 43.6, 46.5); arrow(37.4, 41.5, 43.6, 42.5); arrow(37.4, 28, 43.6, 38.5)
    arrow(58.4, 45, 64.6, 54); arrow(72, 49.6, 72, 46.4); arrow(79.4, 41.5, 83.6, 41.5)
    arrow(72, 36.6, 72, 32.4); arrow(91.5, 36.6, 91.5, 32.4); arrow(79.4, 27, 83.6, 27)
    for x, y, t in ((42.6, 51.4, "raw\nMAC"), (40.5, 43.2, "shape"), (43.0, 30.8, "magni-\ntude"), (59.6, 51.2, "calibrated\nMAC")):
        ax.text(x, y, t, fontsize=5.2, color=C["ink2"], ha="center", va="center", linespacing=0.95)
    for kind, lab, x in (("computed", "computed", 2), ("measured", "measured input", 18), ("model", "model step", 38),
                         ("obs", "data role / outcome", 55)):
        ax.add_patch(FancyBboxPatch((x, 60.2), 2.2, 1.4, boxstyle="round,pad=0.1", fc=style[kind][0], ec=style[kind][1], lw=0.6))
        ax.text(x + 3, 60.9, lab, fontsize=6, va="center")
    data = pd.DataFrame(dict(element=["TD-DFT", "HPLC shape", "extract MAC", "calibration", "cell optics", "BioSNICAR",
                                      "observation operator", "host interface"],
                             status=["computed", "measured (calibration target)", "measured (calibration target)", "fitted",
                                     "model + measured inputs", "model", "model (scalar k)", "offline reference SEB"]))
    save(fig, "Fig1_workflow", data,
         "Workflow from quantum-chemical pigment spectra to glacier radiation. Orange: computed; green: measured "
         "inputs (both are CALIBRATION TARGETS, not validation data); blue: model steps; grey: data roles and "
         "outcomes. Only the calibrated spectral shape reaches the glacier model (the TD-DFT oscillator-strength "
         "magnitude is absorbed by the fitted scale f together with the phenol-equivalent assay conversion). The "
         "observation operator uses a scalar HCRF/albedo factor k whose prior overlaps the S6 test plots (D3). The "
         "host interface is an offline reference point SEB, not E3SM/MAR.", "schematic (no data)")


# =========================================================================== Figure 2
def load_cals():
    import tddft_calibration as TC  # noqa: F401
    cals = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "phase4/results/cache/calibrations/*.pkl"))):
        c = pickle.load(open(f, "rb"))
        kind = c.diagnostics["residual_model"]["kind"]
        func = "CAM-B3LYP" if "CAM" in c.spec.source else "B3LYP"
        cals[(func, kind)] = (c, os.path.relpath(f, ROOT))
    return cals


def fig2():
    import empirical_data as ED
    import tddft_calibration as TC
    cals = load_cals()
    wl = np.arange(250.0, 750.5, 1.0)
    st_b3 = pd.read_csv(src("phase1/results/level2/level2_B3LYP_states.csv"))
    st_cam = pd.read_csv(src("phase1/results/legacy_2026-10-09/level2_cam/level2_CAM-B3LYP_states.csv"))
    wl_h, S, sS = ED.chromophore_hplc_shape()
    wl_m, E, sE = ED.phenolic_extract_mac()
    src("data/empirical/williamson2020_phenolics_hplc_abs_spec.csv"); src("data/empirical/williamson2020_phenolics_mac_with_error.csv")
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.2))
    rows = []
    for ax, (func, st, col, lab) in zip(axs[0], (("B3LYP", st_b3, C["tddft_D"], "a"), ("CAM-B3LYP", st_cam, C["cam"], "b"))):
        c, cpath = cals[(func, "ar1")]
        USED[cpath] = hashlib.sha256(open(os.path.join(ROOT, cpath), "rb").read()).hexdigest()
        raw = TC.perturbed_mac(c.spec, 0.0, 1.0, 0.3)
        mh = raw(wl_h); shape = mh / np.trapezoid(mh, wl_h)
        ax.fill_between(wl_h, S - 2 * sS, S + 2 * sS, color=C["measured"], alpha=0.25, lw=0, label="HPLC chromophore ±2 SD")
        ax.plot(wl_h, S, color=C["measured"], lw=1.5)
        ax.plot(wl_h, shape, color=col, lw=1.5, label=f"raw TD-{func}\n(FWHM 0.3 eV, dE = 0)")
        ax2 = ax.twinx(); ax2.grid(False); ax2.spines["right"].set_visible(True)
        lam = st["Wavelength_nm"].to_numpy(); f = st["Oscillator_Strength"].to_numpy()
        ax2.vlines(lam, 0, f, color=col, lw=0.8, alpha=0.8); ax2.set_ylabel("oscillator strength (sticks)", color=C["ink2"])
        ax2.set_ylim(0, max(0.7, f.max() * 1.1))
        dE = c.mean("dE")
        ax.axvspan(250, 350, color=C["grid"], alpha=0.6, lw=0)
        ax.text(300, -0.13, "grey: held at\n350 nm in model", transform=ax.get_xaxis_transform(), fontsize=5.2, color=C["ink2"], ha="center", va="top")
        ax.set_xlim(250, 650); ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("unit-area absorbance (nm⁻¹)")
        ax.set_title(f"{func}: raw spectrum vs HPLC shape (calibration target)\n"
                     f"calibrated shift dE = {dE:+.3f} ± {c.sd('dE'):.3f} eV", fontsize=7)
        ax.legend(loc="upper right", fontsize=5.8); panel(ax, lab)
        for x, y in zip(wl_h, shape):
            rows.append(dict(panel=lab, series=f"raw {func} unit-area", wavelength_nm=x, value=y))
        for x, y in zip(lam, f):
            rows.append(dict(panel=lab, series=f"{func} stick f", wavelength_nm=x, value=y))
    for x, y, s in zip(wl_h, S, sS):
        rows.append(dict(panel="a,b", series="HPLC mean shape", wavelength_nm=x, value=y, sd=s))
    ax = axs[1, 0]
    rng = np.random.default_rng(0)
    for func, col in (("B3LYP", C["tddft_D"]), ("CAM-B3LYP", C["cam"])):
        c, _ = cals[(func, "ar1")]
        idx = rng.choice(c.samples.shape[0], 200, replace=False)
        D = np.array([TC.complexed_mac(c.spec, *c.samples[i, [0, 2, 1, 3]][[0, 1, 2, 3]], )(wl) if False else
                      TC.complexed_mac(c.spec, c.samples[i, 0], c.samples[i, 2], c.samples[i, 1], c.samples[i, 3])(wl)
                      for i in idx])
        Cc = np.array([TC.perturbed_mac(c.spec, c.samples[i, 0], c.samples[i, 2], c.samples[i, 1])(wl) for i in idx])
        q = np.quantile(D, [0.025, 0.5, 0.975], axis=0)
        ax.fill_between(wl, q[0], q[2], color=col, alpha=0.18, lw=0)
        ax.plot(wl, q[1], color=col, lw=1.5, label=f"{func} D (median, 95 % of draws)")
        ax.plot(wl, np.median(Cc, axis=0), color=col, lw=1.2, ls="--", label=f"{func} C (no Fe term)")
        for x, a, b, m in zip(wl[::5], q[0][::5], q[2][::5], q[1][::5]):
            rows.append(dict(panel="c", series=f"{func} D", wavelength_nm=x, value=m, lo=a, hi=b))
    ax.fill_between(wl_m, E - 2 * sE, E + 2 * sE, color=C["measured"], alpha=0.25, lw=0)
    ax.plot(wl_m, E, color=C["measured"], lw=1.5, label="extract MAC ±2 SE (target)")
    ax.set_yscale("log"); ax.set_xlim(260, 750); ax.set_ylim(1e2, 2e6)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("MAC (m² kg⁻¹ phenol eq., log scale)")
    ax.set_title("Calibrated MAC vs extract (both fitted to it; not validation)", fontsize=7)
    ax.legend(fontsize=5.3, loc="lower left", frameon=True, facecolor="white", edgecolor="none", framealpha=0.92); panel(ax, "c")
    ax = axs[1, 1]
    for func, col in (("B3LYP", C["tddft_D"]), ("CAM-B3LYP", C["cam"])):
        c, _ = cals[(func, "ar1")]
        ax.hist(c.samples[:, 0], bins=50, color=col, alpha=0.75, label=f"{func}: dE posterior")
    ax.axvline(0, color=C["ink2"], lw=0.6)
    ax.set_xlabel("calibrated band shift dE (eV)"); ax.set_ylabel("posterior samples")
    ax.set_title("The calibration absorbs a 0.8 eV functional error", fontsize=7)
    ax.legend(fontsize=6); panel(ax, "d")
    fig.tight_layout()
    save(fig, "Fig2_spectra", pd.DataFrame(rows),
         "Raw and calibrated pigment spectra. (a, b) Raw TD-B3LYP and TD-CAM-B3LYP spectra (sticks: oscillator "
         "strengths, right axis; curve: unit-area broadened absorbance, FWHM 0.3 eV, no shift) against the "
         "isolated-chromophore HPLC shape (mean ± 2 SD over HPLC peaks 2–4, Williamson et al. 2020). Grey: < 350 nm, "
         "where the glacier model holds the MAC at its 350 nm value. (c) Calibrated MACs: median and 95 % range of "
         "200 joint posterior draws (AR(1) residual model), tier D solid, tier C dashed, against the measured whole-"
         "extract MAC ± 2 SE (53 S6 samples). Log axis: the MAC spans 4 orders of magnitude. Both the HPLC shape and "
         "the extract MAC are calibration targets; agreement with them is NOT validation. (d) Posterior of the band "
         "shift dE: B3LYP +0.04 eV, CAM-B3LYP −0.79 eV.", "calibration-data comparison (not independent validation)",
         dict(calibration="AR(1) cut posterior, cached by content fingerprint", draws=200))


# =========================================================================== Figure 3
def fig3():
    p = "records/molecular_contribution/uncertainty/posterior_contrasts.csv"
    if not os.path.isfile(os.path.join(ROOT, p)):
        return blocked("Fig3_substitution", "posterior contrasts not computed")
    P = pd.read_csv(src(p))
    st = json.load(open(src("records/molecular_contribution/uncertainty/settings.json")))
    pairs = [("T-B3-D - T-CAM-D", "B3LYP-D − CAM-D (independent posteriors)", C["cam"]),
             ("T-B3-D - T-MEAS", "B3LYP-D − measured MAC", C["measured"]),
             ("T-CAM-D - T-MEAS", "CAM-D − measured MAC", C["measured"]),
             ("T-B3-D - T-B3-C", "B3LYP-D − B3LYP-C (Fe/visible term; paired)", C["tddft_D"])]
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True)
    P = P.sort_values(["log_b", "r_um", "dust_ppb"])
    states = P[["log_b", "r_um", "dust_ppb"]].drop_duplicates().reset_index(drop=True)
    states["x"] = np.arange(len(states))
    P = P.merge(states, on=["log_b", "r_um", "dust_ppb"])
    rows = []
    for ax, q, thr, unit in ((axs[0], "bba", 0.01, "Δ broadband albedo"), (axs[1], "absorbed_algal_W", 10.0, "Δ absorbed SW (W m⁻²)")):
        for k, (pair, lab, col) in enumerate(pairs):
            d = P[(P.pair == pair) & (P.quantity == q)]
            off = (k - 1.5) * 0.18
            ls = "--" if "CAM-D − measured" in lab else "-"
            ax.vlines(d.x + off, d.q025, d.q975, color=col, lw=1.4, linestyles=ls, label=f"{lab}: 95 % posterior")
            ax.plot(d.x + off, d.plug_in, "o", ms=3.2, mfc="white", mec=col, mew=0.9)
            for r in d.itertuples():
                rows.append(dict(quantity=q, pair=pair, log_b=r.log_b, r_um=r.r_um, dust_ppb=r.dust_ppb, plug_in=r.plug_in,
                                 q025=r.q025, q975=r.q975, p_exceeds_threshold=r.p_exceeds_threshold,
                                 measured_mac_uncertainty=r.measured_mac_uncertainty))
        ax.axhspan(-thr, thr, color=C["grid"], alpha=0.7, lw=0)
        ax.axhline(0, color=C["ink2"], lw=0.6)
        ax.set_ylabel(unit)
        ax.text(-0.6, thr, f"practical threshold ±{thr:g}", fontsize=5.8, ha="left", va="bottom", color=C["ink2"])
    axs[0].plot([], [], "o", ms=3.2, mfc="white", mec=C["ink2"], label="plug-in (posterior-mean spectra)")
    axs[0].legend(fontsize=5.6, ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.01), frameon=False)
    axs[1].set_xticks(states.x)
    axs[1].set_xticklabels([f"{b:g}|{r / 1000:g}|{d / 1e5:g}" for b, r, d in states[["log_b", "r_um", "dust_ppb"]].to_numpy()],
                           rotation=90, fontsize=5.5)
    axs[1].set_xlabel("state: log₁₀ cells mL⁻¹ | bubble radius (mm) | dust (10⁵ ppb)")
    panel(axs[0], "a"); panel(axs[1], "b")
    fig.tight_layout()
    meas = P[P.pair.str.contains("MEAS")].measured_mac_uncertainty.iloc[0]
    save(fig, "Fig3_substitution", pd.DataFrame(rows),
         f"Controlled pigment-level substitution (Experiment A). Only the phenolic absorption spectrum changes; cells, "
         f"packaging, photosynthetic pigments, column, bubbly ice, dust and illumination (SZA 50°, clear sky) are "
         f"identical. Bars: 95 % posterior intervals of the signed contrast (TD-DFT: {st['n_draws_per_TD_treatment']} "
         f"joint calibration draws per treatment; B3LYP vs CAM-B3LYP: independent posterior draws, all pairs; C vs D: "
         f"paired by posterior sample; measured-MAC uncertainty: {meas}). Open circles: plug-in contrasts from "
         f"posterior-mean spectra. Grey band: practical-importance thresholds (0.01 albedo, 10 W m⁻²), which are not "
         f"measurement uncertainties. The Fe/visible-absorber term (orange, D − C) is an empirical ablation (the Fe "
         f"increment is measured, not computed). 24 states.", "model experiment (no observations)",
         dict(n_draws=st["n_draws_per_TD_treatment"], thresholds=st["thresholds"]))


# =========================================================================== Figure 4
def fig4():
    R = pd.read_csv(src("records/forward_diagnostics/per_plot.csv"))
    O = pd.read_csv(src("records/forward_diagnostics/observations.csv")).set_index("sample")
    RS = pd.read_csv(src("records/forward_diagnostics/residual_spectra_5nm.csv"))
    SB = pd.read_csv(src("records/molecular_contribution/albedo_slope_day_bootstrap.csv"))
    UF = pd.read_csv(src("records/molecular_contribution/unexplained_fraction_day_bootstrap.csv"))
    setj = json.load(open(src("records/forward_diagnostics/settings.json")))
    ref = R[R.variant == "reference"]
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 3)
    rows = []
    models = ["no_algae", "tierA_empirical", "tddft_D"]
    for k, m in enumerate(models):
        ax = fig.add_subplot(gs[0, k])
        d = ref[ref.optics == m]
        zero = d.cells <= 0
        ax.plot([0, 0.8], [0, 0.8], color=C["ink2"], lw=0.6)
        ax.plot(d.alb_obs[~zero], d.alb_mod[~zero], "o", ms=3.5, mfc=COL[m], mec="white", mew=0.5, label="counted plots")
        ax.plot(d.alb_obs[zero], d.alb_mod[zero], "s", ms=4, mfc="white", mec=C["ink"], mew=0.9, label="zero-count plots")
        ax.set_xlim(0, 0.8); ax.set_ylim(0, 0.8); ax.set_aspect("equal")
        ax.set_xlabel("observed broadband albedo"); ax.set_ylabel("predicted (reference state)" if k == 0 else "")
        ax.set_title(f"{NAMES[m]}\nbias {d.err_alb.mean():+.3f}, MAE {d.err_alb.abs().mean():.3f}", fontsize=7)
        if k == 0:
            ax.legend(fontsize=5.8, loc="upper left")
        panel(ax, "abc"[k])
        for r in d.itertuples():
            rows.append(dict(panel="abc"[k], optics=m, sample=r.sample, cells=r.cells, observed=r.alb_obs, predicted=r.alb_mod))
    ax = fig.add_subplot(gs[1, 0:2])
    for m in ("no_algae", "tierA_empirical", "tddft_D", "measured_mac_C"):
        if f"{m}_median" not in RS:
            continue
        ax.fill_between(RS.wavelength_nm, RS[f"{m}_q25"], RS[f"{m}_q75"], color=COL[m], alpha=0.15, lw=0)
        ax.plot(RS.wavelength_nm, RS[f"{m}_median"], color=COL[m], lw=1.4, ls=LS.get(m, "-"), label=NAMES[m])
        for x, y in zip(RS.wavelength_nm, RS[f"{m}_median"]):
            rows.append(dict(panel="d", optics=m, wavelength_nm=x, residual_median=y))
    ax.axhline(0, color=C["ink2"], lw=0.6)
    for b, w in (("B2", 490), ("B4", 665), ("B8", 833)):
        ax.axvline(w, color=C["grid"], lw=1.0); ax.text(w + 8, 0.54, b, fontsize=6, color=C["ink2"])
    ax.set_xlim(350, 1300); ax.set_ylim(-0.3, 0.6)
    ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("model − observed albedo\n(median, IQR over 46 plots)")
    ax.set_title("Spectral residual at known abundance (reference nuisance state)", fontsize=7)
    ax.legend(fontsize=5.8, ncol=2); panel(ax, "d")
    ax = fig.add_subplot(gs[1, 2])
    x = np.arange(3)
    so = SB.set_index("band").loc[["B2", "B4", "B8"]]
    ax.bar(x - 0.18, so.slope_obs, 0.34, color=C["ink2"], label="observed")
    ax.errorbar(x - 0.18, so.slope_obs, yerr=[so.slope_obs - so.obs_lo, so.obs_hi - so.slope_obs], fmt="none", ecolor=C["ink"], lw=0.8, capsize=2)
    mod = so.slope_obs - so.unexplained
    ax.bar(x + 0.18, mod, 0.34, color=C["tddft_D"], label="TD-DFT D model")
    ax.set_xticks(x); ax.set_xticklabels(["B2 490", "B4 665", "B8 833"]); ax.set_ylabel("Δ albedo per decade of cells")
    uf = UF[(UF.optics == "tddft_D") & (UF.band == "broadband")].iloc[0]
    ax.set_ylim(min(so.obs_lo.min(), mod.min()) * 1.12, 0)
    ax.set_title(f"Slope vs abundance (41 plots)\nunexplained broadband\nfraction {uf.unexplained_fraction:.2f} "
                 f"({uf.frac_q025:.2f}–{uf.frac_q975:.2f})", fontsize=6.5)
    ax.legend(fontsize=5.8, loc="lower right"); panel(ax, "e")
    for b in ("B2", "B4", "B8"):
        rows.append(dict(panel="e", band=b, slope_obs=so.loc[b, "slope_obs"], obs_lo=so.loc[b, "obs_lo"], obs_hi=so.loc[b, "obs_hi"],
                         slope_model=mod.loc[b], unexplained=so.loc[b, "unexplained"], unexpl_lo=so.loc[b, "unexpl_lo"],
                         unexpl_hi=so.loc[b, "unexpl_hi"]))
    fig.tight_layout()
    rs = setj["reference_state"]
    save(fig, "Fig4_forward_diagnosis", pd.DataFrame(rows),
         f"Forward-error diagnosis at known abundance, S6 2017 (development data, previously used for σ, radius prior "
         f"and the k prior). Reference nuisance state fixed a priori: bubble radius {rs['r_um']:.0f} µm, dust "
         f"{rs['dust_ppb']:.3g} ppb, f_n {rs['f_n']:.2f}, crust 2 cm at 450 kg m⁻³. (a–c) Predicted vs observed "
         f"hemispherical broadband albedo (46 plots; 1 excluded for > 5 % missing HCRF; 5 zero-count plots as squares). "
         f"(d) Spectral residual (median and IQR). (e) Albedo change per decade of counted cells: observed (95 % "
         f"day-block bootstrap, 8 days) vs the TD-DFT D model at the reference state (posterior-mean optics, so the model slope carries no interval; the unexplained-fraction interval in the title is the day-block bootstrap). Panel d y-range is limited to −0.3…0.6; no data in 350–1300 nm fall outside it. The red/NIR shortfall is "
         f"darkening associated with abundance that the reference model does not reproduce; it is not attributed to a "
         f"specific cause here (crust water, ice structure, co-located impurities and plot heterogeneity remain "
         f"candidates). Tested geometry adjustments (spectrally flat HCRF/albedo, 0.83–0.88) do not explain the full "
         f"reference-model bias.", "development diagnostics (not independent validation)",
         dict(reference_state=rs))


# =========================================================================== Figure 5
def fig5():
    H1 = pd.read_csv(src("records/h1_albedo/h1_contrasts.csv"))
    H24 = pd.read_csv(src("records/h2h4/h2h4_contrasts.csv"))
    H5 = pd.read_csv(src("records/h5/h5_summary_test.csv")).query("split == 'final_test'")
    src("records/h1_albedo/h1_summary.csv"); src("records/h2h4/h2h4_summary.csv")
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1, 1.15, 0.8]))
    rows = []
    ax = axs[0]
    lab = {"MAE(tierA_empirical) - MAE(tddft_D)": "M1 Tier A − M3", "MAE(measured_mac_C) - MAE(tddft_D)": "M1b measured − M3",
           "MAE(no_algae) - MAE(tddft_D)": "M0 no algae − M3"}
    for i, r in enumerate(H1.itertuples()):
        cm = COL[r.contrast.split("(")[1].split(")")[0]]
        ax.errorbar(r.mean, i, xerr=[[r.mean - r.lo], [r.hi - r.mean]], fmt="o", ms=4, color=cm, ecolor=cm, capsize=2, lw=1.2)
        rows.append(dict(test="H1", contrast=r.contrast, estimate=r.mean, lo=r.lo, hi=r.hi, n=r.n, blocks=r.n_blocks, units="albedo MAE"))
    ax.set_yticks(range(len(H1))); ax.set_yticklabels([lab.get(c, c) for c in H1.contrast], fontsize=6)
    ax.axvline(0, color=C["ink2"], lw=0.6); ax.axvline(0.01, color=C["ink2"], lw=0.8, ls="--")
    ax.set_xlabel("MAE difference (albedo)\n> 0: M3 better")
    ax.set_title(f"H1 plot albedo (S6 2017)\nn = {int(H1.n.iloc[0])} plots, {int(H1.n_blocks.iloc[0])} day-blocks", fontsize=7)
    panel(ax, "a")
    ax = axs[1]
    P = H24[H24.set == "primary"].reset_index(drop=True)
    labs = []
    for i, r in enumerate(P.itertuples()):
        thr = 10.0 if r.hypothesis == "H2" else 0.01
        scale = 1.0 if r.hypothesis == "H2" else 1000.0
        cmp = r.contrast.split(": ")[1].split(" ")[0]
        col = {"M0": C["none"], "M1": C["tierA"], "M1b": C["measured"], "M2": C["m2"]}.get(cmp, C["ink2"])
        ax.errorbar(r.mean * scale, i, xerr=[[(r.mean - r.lo) * scale], [(r.hi - r.mean) * scale]], fmt="o", ms=4, color=col, capsize=2, lw=1.2)
        labs.append(f"{r.hypothesis}: {r.contrast.split(': ')[1]}" + (" (×10⁻³ albedo)" if r.hypothesis == "H4" else " (W m⁻²)"))
        rows.append(dict(test=r.hypothesis, contrast=r.contrast, estimate=r.mean, lo=r.lo, hi=r.hi, n=r.n, blocks=r.n_blocks,
                         units="W m-2 MAE" if r.hypothesis == "H2" else "albedo MAE"))
    ax.set_yticks(range(len(P))); ax.set_yticklabels(labs, fontsize=5.6)
    ax.axvline(0, color=C["ink2"], lw=0.6); ax.axvline(10, color=C["ink2"], lw=0.8, ls="--")
    ax.set_xlabel("MAE difference, > 0: M3 better\n(dashed: 10 W m⁻² or 10×10⁻³ albedo)")
    ax.set_title(f"H2/H4 KAN_L/KAN_M overpass windows\nn = {int(P.n.iloc[0])} station-days, {int(P.n_blocks.iloc[0])} station-years", fontsize=7)
    panel(ax, "b")
    ax = axs[2]
    H5 = H5.set_index("model")
    order = ["no_algae", "tierA_empirical", "measured_mac_C", "tddft_D", "tddft_C"]
    for i, m in enumerate(order):
        ax.barh(i, H5.loc[m, "mae"], color=COL[m], height=0.6, hatch="////" if m == "tddft_C" else None, edgecolor="white", lw=0)
        ax.text(H5.loc[m, "mae"] + 0.002, i, f"{H5.loc[m, 'mae']:.3f}", va="center", fontsize=5.8)
        rows.append(dict(test="H5", model=m, mae=H5.loc[m, "mae"], n=int(H5.loc[m, "n"]), units="albedo MAE"))
    ax.set_yticks(range(len(order))); ax.set_yticklabels([NAMES[m] for m in order], fontsize=5.8)
    ax.set_xlabel("albedo MAE"); ax.set_xlim(0, 0.16)
    ax.set_title("H5 (PROMBIO)\n2 test days:\ninsufficient", fontsize=7)
    panel(ax, "c")
    fig.tight_layout()
    save(fig, "Fig5_benchmark", pd.DataFrame(rows),
         "Glacier-model benchmark against all prespecified baselines (frozen protocol, amendments A1, A2/A2.1). "
         "(a) H1: plot-scale broadband albedo, S6 2017 (47 plots, 9 sampling-day blocks; leave-one-day-out empirical "
         "Bayes; development-reused sites). (b) H2 (absorbed SW, W m⁻²) and H4 (albedo, ×10⁻³) at KAN_L/KAN_M over "
         "±1 h of the Sentinel-2 overpass (126 station-days, 11 station-year blocks; 2019 reported separately). M2 is "
         "a simple four-band irradiance-weighted conversion (protocol fallback), NOT a published method. (c) H5: station "
         "albedo from PROMBIO biology under the radiometer (2 test station-days; no interval). Error bars: 95 % "
         "cluster-bootstrap intervals over the stated blocks. Dashed: prespecified minimum meaningful improvements. "
         "Point colour = the comparator model (M0 grey, M1 blue, M1b green, M2 pink). H1 intervals for M1 and M1b are narrower than the markers (≤ ±0.0004). Hatched bar: tier C (no Fe term). Every M3 superiority hypothesis is NOT SUPPORTED. Deviations: posterior-mean optics (D1); k-prior overlap "
         "with S6 plots (D3, H1 does not use k).", "independent tests per protocol (H1 on previously seen sites)")


# =========================================================================== Figure 6
def fig6():
    cand = glob.glob(os.path.join(ROOT, "phase4/results/scene_product/*/*_provenance.json"))
    if not cand:
        return blocked("Fig6_satellite", "no completed, validated scene product exists (scene_product.py never ran to "
                                         "completion; repeated container reboots). A scene figure from legacy maps "
                                         "(phase4/maps4.py, pre-repair) would be stale and is not produced.")


# =========================================================================== Supplementary
def figS1():
    import re
    rows = []
    for job in ("L1_FULL", "L2_CAM_TDA15"):
        log = open(src(f"phase1/results/v2/{job}/job.log"), errors="replace").read()
        for m in re.finditer(r"residual target ([0-9.e-]+): (\d+)/(\d+) roots converged after <= (\d+) cycles \((\d+) s this chunk\)", log):
            rows.append(dict(job=job, target=float(m.group(1)), converged=int(m.group(2)), roots=int(m.group(3)),
                             cycles=int(m.group(4)), chunk_s=int(m.group(5))))
    D = pd.DataFrame(rows)
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for job, col in (("L1_FULL", C["tierA"]), ("L2_CAM_TDA15", C["cam"])):
        d = D[D.job == job]
        axs[0].plot(d.cycles, d.converged / d.roots, "o-", color=col, ms=3.5, label=f"{job} (residual 0.01 stage)")
        axs[1].plot(d.cycles, d.chunk_s / 60, "o-", color=col, ms=3.5, label=job)
    axs[0].set_xlabel("cumulative Davidson cycles"); axs[0].set_ylabel("fraction of roots converged"); axs[0].set_ylim(0, 1)
    axs[0].legend(fontsize=6); panel(axs[0], "a")
    axs[1].set_xlabel("cumulative Davidson cycles"); axs[1].set_ylabel("wall time per 2-cycle chunk (min)")
    axs[1].axhline(32, color=C["ink2"], ls="--", lw=0.8); axs[1].text(1, 33, "shortest observed uptime (32 min)", fontsize=5.8, color=C["ink2"])
    axs[1].legend(fontsize=6); panel(axs[1], "b")
    fig.tight_layout()
    save(fig, "FigS1_chemistry_progress", D,
         "Chemistry convergence progress from durable checkpoints (job logs). Only the first stage (residual 0.01) "
         "has been reached by either job; neither has completed a stage. L1_FULL: 18 cycles, 24/30 roots; "
         "L2_CAM_TDA15: 2 cycles, 6/15 roots (one 90-min chunk). No checkpoint advanced across the 19:52 and 20:57 UTC "
         "relaunches. Matched TD vs TDA and root-count checks are not yet available (jobs queued).",
         "progress record (no scientific result)")


def figS2():
    cals = load_cals()
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
    rows = []
    for (func, kind), (c, cpath) in sorted(cals.items()):
        USED[cpath] = hashlib.sha256(open(os.path.join(ROOT, cpath), "rb").read()).hexdigest()
        col = C["cam"] if func == "CAM-B3LYP" else C["tddft_D"]
        mk = "o" if kind == "ar1" else "^"
        s = c.samples[::8]
        axs[0].plot(s[:, 0], s[:, 1], mk, ms=1.6, color=col, alpha=0.5, label=f"{func} {kind}")
        axs[1].plot(np.log(s[:, 2]), s[:, 3], mk, ms=1.6, color=col, alpha=0.5, label=f"{func} {kind}")
        d = c.diagnostics
        rows.append(dict(functional=func, residuals=kind, rho_shape=d["residual_model"].get("rho_shape_mean"),
                         mac_r2_log_400_700=d.get("mac_r2_log_400_700"), fe_share_vis=d.get("fe_share_of_absorption_400_700"),
                         at_bound=json.dumps(d.get("at_bound")), corr=json.dumps(c.summary()["corr_dE_w_f_phi"])))
    axs[0].set_xlabel("dE (eV)"); axs[0].set_ylabel("w (eV)"); axs[0].legend(fontsize=5.6, markerscale=3); panel(axs[0], "a")
    axs[1].set_xlabel("ln f"); axs[1].set_ylabel("φ (Fe fraction, effective)"); panel(axs[1], "b")
    T = pd.DataFrame(rows)
    ax = axs[2]
    x = np.arange(len(T))
    ax.bar(x, T.mac_r2_log_400_700, color=[C["cam"] if f == "CAM-B3LYP" else C["tddft_D"] for f in T.functional],
           hatch=["////" if k == "iid" else None for k in T.residuals], edgecolor="white", lw=0)
    ax.set_xticks(x); ax.set_xticklabels([f"{f}\n{k}" for f, k in zip(T.functional, T.residuals)], fontsize=5.6)
    ax.axhline(0, color=C["ink2"], lw=0.6); ax.set_ylabel("log-R² of extract MAC, 400–700 nm")
    ax.set_title("visible fit depends on residual model", fontsize=7); panel(ax, "c")
    fig.tight_layout()
    save(fig, "FigS2_calibration", T,
         "Calibration posterior dependence and adequacy. (a) Shift–width and (b) scale–Fe-fraction joint posteriors "
         "(every 8th of 4000 samples; circles AR(1), triangles independent residuals). φ is an effective mixture "
         "coefficient, not a measured complexed fraction. (c) Visible-region (400–700 nm) log-R² of the calibrated "
         "tier-D MAC against the extract it was fitted to: −0.07 (B3LYP, AR(1); shape-residual ρ ≈ 0.94), 0.67 "
         "(B3LYP, iid), 0.19 (CAM-B3LYP, AR(1)). The visible fit, which carries the glacier-relevant absorption, is "
         "poor and residual-model dependent; a refined quadrature does not make the residual model adequate.",
         "calibration diagnostics")


def figS3():
    S = pd.read_csv(src("records/ml_comparison_physics/summary.csv"))
    H1 = pd.read_csv(src("records/h1_albedo/h1_summary.csv")).set_index("model")
    EN = pd.read_csv(src("records/forward_diagnostics/ensemble_summary.csv")).set_index("optics")
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8))
    ax = axs[0]
    ms = ["no_algae", "tierA_empirical", "measured_mac_C", "tddft_D", "tddft_C"]
    for k, m in enumerate(ms):
        ax.plot(H1.loc[m, "width90"], H1.loc[m, "coverage90"], "o", ms=5, color=COL[m], mfc=COL[m] if m != "tddft_C" else "white",
                alpha=0.85, label=NAMES[m])
        ax.plot(EN.loc[m, "width90"], EN.loc[m, "coverage90"], "s", ms=4.5, color=COL[m], mfc="white")
    ax.axhline(0.9, color=C["ink2"], ls="--", lw=0.7)
    ax.set_xlabel("mean 90 % interval width (albedo)"); ax.set_ylabel("empirical 90 % coverage")
    ax.set_title("Coverage vs sharpness: H1 (●, discrepancy SD by EB)\nvs nuisance prior only (□)", fontsize=7)
    ax.set_ylim(0, 1); ax.legend(fontsize=5.6, loc="center", title="optics (colour)", title_fontsize=5.6); panel(ax, "a")
    ax = axs[1]
    rows = []
    for fold, off in (("primary", -0.2), ("secondary", 0.2)):
        d = S[S.fold == fold].reset_index(drop=True)
        cols = [C["ink2"] if m.startswith("ml_") else COL.get(m.replace("physics_", ""), C["ink2"]) for m in d.model]
        ax.barh(np.arange(len(d)) + off, d.rmse_pooled, height=0.38, color=cols, alpha=0.9 if fold == "primary" else 0.5)
        for r in d.itertuples():
            rows.append(dict(fold=fold, model=r.model, rmse_pooled=r.rmse_pooled))
    ax.set_yticks(np.arange(len(d))); ax.set_yticklabels([("physics: " + NAMES.get(m.replace("physics_", ""), m)) if m.startswith("physics_") else m.replace("ml_", "ML: ").replace("_", " ") for m in d.model], fontsize=5.6)
    ax.set_xlabel("pooled 4-band HCRF RMSE"); ax.set_title("Forward: physics vs ML, same rows\n(solid primary fold, faded secondary)", fontsize=7)
    panel(ax, "b")
    fig.tight_layout()
    D = pd.concat([pd.DataFrame(rows), H1[["coverage90", "width90"]].reset_index().assign(kind="H1"),
                   EN[["coverage90", "width90"]].reset_index().rename(columns={"optics": "model"}).assign(kind="prior_ensemble")])
    save(fig, "FigS3_coverage_ml", D,
         "(a) Predictive coverage vs sharpness for broadband albedo at S6 2017: H1 predictive distributions (filled; "
         "discrepancy SD chosen by empirical Bayes on training days) reach about 0.9 coverage only with wide intervals, "
         "while the nuisance-prior-only ensemble (open squares) covers 22–24 % (biased and too narrow). (b) Forward 4-band "
         "HCRF RMSE on identical positive-count rows: physics point predictions (colour) vs ML baselines (grey) from "
         "records/ml_comparison (retrospective, development data). Different endpoints from (a). In (a) Measured MAC "
         "(M1b) and TD-DFT D (M3) coincide (coverage 0.894, width 0.366) and their markers overlap; Tier A, M1b and M3 "
         "squares also overlap. Exact values are in the data CSV.",
         "development diagnostics")


def figS4():
    S = pd.read_csv(src("records/forward_diagnostics/summary.csv"))
    d = S[S.optics == "tddft_D"].reset_index(drop=True)
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw=dict(width_ratios=[1.4, 1]))
    ax = axs[0]
    y = np.arange(len(d))
    ax.barh(y - 0.2, d.bias_alb, 0.38, color=C["tddft_D"], label="broadband bias")
    ax.barh(y + 0.2, d.bias_B2, 0.38, color=C["tierA"], label="B2 bias")
    ax.set_yticks(y); ax.set_yticklabels(d.variant, fontsize=5.5); ax.axvline(0, color=C["ink2"], lw=0.6)
    ax.set_xlabel("model − observed albedo (46 plots)"); ax.legend(fontsize=6); ax.invert_yaxis(); panel(ax, "a")
    ax = axs[1]
    mc = os.path.join(ROOT, "records/molecular_contribution/uncertainty/mc_stability_24_vs_all.csv")
    rows = d[["variant", "bias_alb", "bias_B2"]].assign(panel="a")
    if os.path.isfile(mc):
        M = pd.read_csv(src("records/molecular_contribution/uncertainty/mc_stability_24_vs_all.csv"))
        M = M[M.quantity.isin(["bba"])]
        PL = {"T-B3-C - T-CAM-C": "B3LYP-C − CAM-C", "T-B3-D - T-B3-C": "B3LYP D − C", "T-B3-D - T-CAM-D": "B3LYP-D − CAM-D",
              "T-B3-D - T-MEAS": "B3LYP-D − measured*", "T-CAM-D - T-CAM-C": "CAM D − C", "T-CAM-D - T-MEAS": "CAM-D − measured*"}
        y = np.arange(len(M))
        ax.barh(y - 0.19, M.max_d_q025, 0.36, color=C["none"], label="lower (2.5 %) end")
        ax.barh(y + 0.19, M.max_d_q975, 0.36, color=C["ink2"], label="upper (97.5 %) end")
        ax.set_yticks(y); ax.set_yticklabels([PL.get(q, q) for q in M.pair], fontsize=5.6)
        ax.set_xlabel("max shift of interval end over states,\n24 vs 96 draws (broadband albedo)")
        ax.legend(fontsize=5.6, loc="lower right")
        ax.set_title("Monte Carlo stability", fontsize=7)
        rows = pd.concat([rows, M.assign(panel="b")])
    else:
        ax.text(0.5, 0.5, "MC stability not yet computed", ha="center", transform=ax.transAxes); ax.axis("off")
    panel(ax, "b")
    fig.tight_layout()
    save(fig, "FigS4_sensitivity", rows,
         "(a) One-factor-at-a-time sensitivity of the TD-DFT D forward bias at known abundance (S6 2017, development "
         "data; reference state as Fig. 4). Ranges: radius and dust at prior 10/90 %, crust density Cooper et al. "
         "(2018) range, crust depth 1/5 cm, species fraction 0/1, algae in a 1 mm film, granular ice, SZA ± 5°, fully "
         "diffuse, mid-latitude spectrum. Only absorber amount per column (crust depth, dust) substantially reduces the "
         "visible bias; none removes it. (b) Monte Carlo stability of the posterior contrasts: largest shift over the 24 states of each 95 % interval end "
         "between 24 and 96 draws per treatment. *The 96-draw measured-MAC contrasts also add the measured-MAC "
         "uncertainty, so their shift combines Monte Carlo error with the added uncertainty source.",
         "development diagnostics / numerical check")


def figS5():
    S = pd.read_csv(src("records/host_coupling/states_and_optics.csv"))
    R = pd.read_csv(src("records/host_coupling/coupling_runs.csv"))
    num = json.load(open(src("records/host_coupling/numerics.json")))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6))
    d = pd.to_datetime(S.day)
    axs[0].plot(d, S.dalpha_tierA_empirical, "o--", color=C["tierA"], ms=3.5, lw=0.8, label="Tier A Δα")
    axs[0].plot(d, S.dalpha_tddft_D, "o--", color=C["tddft_D"], ms=3.5, lw=0.8, label="TD-DFT D Δα")
    axs[0].set_ylabel("algal albedo reduction Δα"); axs[0].legend(fontsize=6); axs[0].tick_params(axis="x", labelsize=5.5, rotation=30)
    panel(axs[0], "a")
    x = np.arange(len(R))
    axs[1].bar(x, R.modelled_increment_mwe * 1000, color=[C["tierA"] if "tierA" in r else C["tddft_D"] for r in R.run])
    axs[1].set_xticks(x); axs[1].set_xticklabels([r.replace("modeA_", "Mode A\n").replace("modeB_", "Mode B\n").replace("tierA_empirical", "Tier A").replace("tddft_D", "TD-DFT D") for r in R.run], fontsize=5.8)
    axs[1].set_ylabel("modelled algal melt increment\n(mm w.e., 34 days)")
    panel(axs[1], "b")
    fig.tight_layout()
    save(fig, "FigS5_host_seb", pd.concat([S.assign(panel="a"), R.assign(panel="b")]),
         f"Conditional host response (reference point SEB, KAN_L 2022-07-24 to 08-26, identical forcing). (a) Algal "
         f"albedo reduction from the satellite-retrieved state (8 S2 station-days, interpolated), Tier A vs TD-DFT D at "
         f"identical state; markers are retrieval dates, dashed lines the linear interpolation used by the host run. (b) Modelled melt increments (Modes A and B). These are MODELLED quantities conditional on "
         f"the SEB, forcing and retrieved state; they are NOT validated melt (H3 untested; KAN_L 2022 prerequisites are "
         f"rule-dependent). Numerics: n_sub 1 vs 4 seasonal melt {num['numerics']['melt_n_sub1_mwe']:.4f} vs "
         f"{num['numerics']['melt_n_sub4_mwe']:.4f} m w.e.", "conditional model scenario (not melt validation)")


FIGS = dict(fig1=fig1, fig2=fig2, fig3=fig3, fig4=fig4, fig5=fig5, fig6=fig6, figS1=figS1, figS2=figS2,
            figS3=figS3, figS4=figS4, figS5=figS5)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--only", nargs="+", default=list(FIGS))
    a = p.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    mpath = os.path.join(OUT, "manifest.json")
    old = json.load(open(mpath)) if os.path.isfile(mpath) else {}
    for k in a.only:
        old.pop(k, None)               # drop a stale failure record keyed by function name
        try:
            FIGS[k]()
        except Exception as e:  # noqa: BLE001 - a broken figure is reported, never silently skipped
            MANIFEST[k] = dict(status="FAILED", error=f"{type(e).__name__}: {e}")
            USED.clear()
            print(f"[{k}] FAILED: {e}", flush=True)
    old.update(MANIFEST)
    old["_build"] = dict(build_info(), note="last invocation; each figure entry carries its own build record")
    json.dump(old, open(mpath, "w"), indent=1)
    print(json.dumps({k: v.get("status", "ok") for k, v in MANIFEST.items()}, indent=1))


if __name__ == "__main__":
    main()
