"""
Figures 3A-3C (+ S3 convergence) for the uncertainty / sensitivity analysis.

Shares the Phase 2 style (fonts, tier colours, PDF + 600-dpi PNG). Parameter
scale groups get their own three colours (validated colour-blind-safe set) so
they are never confused with model tiers; every bar also carries a text label.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import figures as F  # noqa: E402  (Phase 2 style + tier colours)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.stats import gaussian_kde  # noqa: E402

GROUP_COLORS = {"molecular": "#e87ba4", "cellular": "#008300", "environmental": "#4a3aa7"}
INK, INK_2 = F.INK, F.INK_2
SHORT = {"dE_ev": r"$\Delta E$", "f_scale": r"$f$ scale", "fwhm_ev": "FWHM",
         "cell_length_um": r"$L$", "cell_diameter_um": r"$d$", "c_internal": r"$c_i$",
         "grain_um": r"$r_{grain}$", "rho_top": r"$\rho$", "conc_cells_ml": "abundance",
         "lmct_eps": r"$\varepsilon_{LMCT}$", "lmct_center_nm": r"$\lambda_{LMCT}$"}
S2_DEFAULT = ["cell_length_um", "cell_diameter_um", "c_internal", "grain_um", "rho_top", "conc_cells_ml"]
OUTPUT_LABELS = {
    "rf_A": r"$RF_A$", "rf_B": r"$RF_B$", "rf_C": r"$RF_C$", "rf_D": r"$RF_D$",
    "eff_C": r"$RF_C$ per $10^4$ cells mL$^{-1}$", "eff_D": r"$RF_D$ per $10^4$ cells mL$^{-1}$",
    "d_CB": r"$RF_C-RF_B$ (packaging)", "d_DC": r"$RF_D-RF_C$ (complexation)",
}


def set_style(usetex=False):
    F.set_style(usetex=usetex)


save = F.save


def _tag(ax, t):
    ax.text(-0.14, 1.02, t, transform=ax.transAxes, fontweight="bold", va="bottom")


# --------------------------------------------------------------------------- #
# Figure 3A                                                                     #
# --------------------------------------------------------------------------- #
def fig3a(mc: pd.DataFrame, demo=False):
    """(a) PDFs of instantaneous forcing, tiers A-D; (b) PDFs of forcing efficiency
    (per 10^4 cells mL^-1). Log-x because abundance is log-uniform; KDE in log space.
    Bars under each panel: median (dot) and 95 % interval (P2.5-P97.5)."""
    fig, axes = plt.subplots(1, 2, figsize=(F.DOUBLE_COL, 3.0), constrained_layout=True)
    for ax, prefix, xlabel, tag in [
        (axes[0], "rf", r"Instantaneous forcing $RF$ (W m$^{-2}$)", "a"),
        (axes[1], "eff", r"Forcing efficiency (W m$^{-2}$ per $10^{4}$ cells mL$^{-1}$)", "b")]:
        ymax = 0.0
        stats_rows = []
        for t, (col, mk, ls, lab) in F.TIERS.items():
            x = mc[f"{prefix}_{t}"].to_numpy()
            x = x[x > 0]
            lx = np.log10(x)
            kde = gaussian_kde(lx)
            grid = np.linspace(lx.min() - 0.3, lx.max() + 0.3, 400)
            d = kde(grid)
            ax.plot(10 ** grid, d, color=col, ls=ls, lw=1.6, label=lab)
            ax.fill_between(10 ** grid, d, color=col, alpha=0.10, lw=0)
            ymax = max(ymax, d.max())
            stats_rows.append((t, col, mk, np.median(x), np.percentile(x, 2.5), np.percentile(x, 97.5)))
        for k, (t, col, mk, med, lo, hi) in enumerate(stats_rows):
            y = -ymax * (0.08 + 0.07 * k)
            ax.plot([lo, hi], [y, y], color=col, lw=2.0, solid_capstyle="butt")
            ax.plot(med, y, marker=mk, color=col, ms=5, mec="white", mew=0.8)
            ax.text(hi * 1.08, y, t, va="center", fontsize=plt.rcParams["legend.fontsize"], color=INK)
        ax.set_xscale("log")
        ax.set_ylim(-ymax * 0.40, ymax * 1.08)
        ax.axhline(0, color=INK_2, lw=0.6)
        ax.set_yticks([t for t in ax.get_yticks() if t >= 0])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"Probability density of $\log_{10}$ value")
        _tag(ax, tag)
    axes[0].legend(loc="upper left")
    if demo:
        fig.suptitle("DEMO input - not Phase 1 results", color="#b00020")
    return fig


# --------------------------------------------------------------------------- #
# Figure 3B                                                                     #
# --------------------------------------------------------------------------- #
def _sobol_bars(ax, tab: pd.DataFrame, space, title):
    meta = {p.name: p for p in space.params}
    t = tab.sort_values("ST", ascending=True)
    y = np.arange(len(t))
    cols = [GROUP_COLORS[meta[n].group] for n in t.index]
    ax.barh(y, t.ST.clip(lower=0), height=0.72, color="none", edgecolor=cols, lw=1.2)
    ax.barh(y, t.S1.clip(lower=0), height=0.72, color=cols, alpha=0.85, edgecolor="none")
    ax.errorbar(t.ST.clip(lower=0), y, xerr=t.ST_conf, fmt="none", ecolor=INK_2, elinewidth=0.8, capsize=1.5)
    ax.set_yticks(y)
    ax.set_yticklabels([meta[n].label for n in t.index])
    ax.set_xlim(0, max(1.0, float((t.ST + t.ST_conf).max()) * 1.05))
    ax.set_xlabel(r"Sobol' index ($S_1$ filled, $S_T$ outline)")
    ax.set_title(title, loc="left")
    ax.grid(axis="y", visible=False)


def fig3b(sobol_tabs: dict, group_tabs: dict, space, demo=False):
    """(a,b) parameter-level S1/ST with 95 % bootstrap CIs for two outputs;
    (c) scale-level (grouped) total-effect indices for several outputs."""
    keys = list(sobol_tabs)[:2]
    fig, axes = plt.subplots(1, 3, figsize=(F.DOUBLE_COL, 3.3), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[1, 1, 0.9]))
    for ax, k, tag in zip(axes[:2], keys, "ab"):
        _sobol_bars(ax, sobol_tabs[k], space, OUTPUT_LABELS.get(k, k))
        _tag(ax, tag)
    for ax in axes[1:2]:
        ax.tick_params(axis="y", labelsize=plt.rcParams["ytick.labelsize"] - 0.5)

    ax = axes[2]
    outs = list(group_tabs)
    groups = list(GROUP_COLORS)
    w = 0.8 / len(groups)
    x = np.arange(len(outs))
    for j, g in enumerate(groups):
        st = [group_tabs[o].ST.get(g, np.nan) for o in outs]
        ci = [group_tabs[o].ST_conf.get(g, np.nan) for o in outs]
        s1 = [group_tabs[o].S1.get(g, np.nan) for o in outs]
        xx = x - 0.4 + w * (j + 0.5)
        ax.bar(xx, st, width=w * 0.9, color="none", edgecolor=GROUP_COLORS[g], lw=1.2)
        ax.bar(xx, np.clip(s1, 0, None), width=w * 0.9, color=GROUP_COLORS[g], alpha=0.85)
        ax.errorbar(xx, st, yerr=ci, fmt="none", ecolor=INK_2, elinewidth=0.8, capsize=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels([("$RF_%s$/cell" % o[-1]) if o.startswith("eff") else OUTPUT_LABELS.get(o, o)
                        for o in outs], rotation=35, ha="right")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel(r"Grouped Sobol' index")
    ax.set_title("By scale", loc="left")
    ax.grid(axis="x", visible=False)
    _tag(ax, "c")
    handles = [Patch(facecolor=c, edgecolor=c, label=g) for g, c in GROUP_COLORS.items()]
    fig.legend(handles=handles, loc="outside lower center", ncol=3)
    if demo:
        fig.suptitle("DEMO input - not Phase 1 results", color="#b00020")
    return fig


# --------------------------------------------------------------------------- #
# Figure 3C                                                                     #
# --------------------------------------------------------------------------- #
def fig3c(tornado: pd.DataFrame, mc: pd.DataFrame, s2: pd.DataFrame | None, space,
          tornado_output: str, surface_output: str, s2_output: str, lo_q=0.05, hi_q=0.95, demo=False,
          s2_params=None):
    """(a) tornado: output change when one parameter moves from its P5 to P95 with all
    others at their medians; (b) Monte Carlo interaction surface: mean output binned over
    cell diameter x pigment concentration;
    (c) second-order Sobol' indices S_ij (pairwise interaction shares) among cell size,
    pigment concentration, ice grain size, density and abundance (full matrix in the CSV)."""
    meta = {p.name: p for p in space.params}
    fig = plt.figure(figsize=(F.DOUBLE_COL, 3.1), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.0])

    # (a) tornado
    ax = fig.add_subplot(gs[0])
    t = tornado.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(t))
    ax.barh(y, t.d_low, color=F.SEQ_BLUE[0], height=0.7, label=rf"P{int(lo_q * 100)}")
    ax.barh(y, t.d_high, color=F.SEQ_BLUE[2], height=0.7, label=rf"P{int(hi_q * 100)}")
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([meta[n].label for n in t.parameter])
    ax.set_xlabel(r"$\Delta$" + OUTPUT_LABELS.get(tornado_output, tornado_output) + r" (W m$^{-2}$)")
    ax.text(0.97, 0.30, f"median\n{t.base.iloc[0]:.1f}" + r" W m$^{-2}$", transform=ax.transAxes,
            ha="right", va="bottom", color=INK_2, fontsize=plt.rcParams["legend.fontsize"])
    ax.set_title("One-at-a-time", loc="left")
    ax.legend(loc="lower right", title="parameter at", fontsize=plt.rcParams["legend.fontsize"] - 0.5,
              title_fontsize=plt.rcParams["legend.fontsize"] - 0.5)
    ax.grid(axis="y", visible=False)
    _tag(ax, "a")

    # (b) binned interaction surface from the Monte Carlo sample
    ax = fig.add_subplot(gs[1])
    xb = np.linspace(5, 15, 9)
    yb = np.logspace(1, np.log10(200), 9)
    z = mc[surface_output].to_numpy()
    H, _, _ = np.histogram2d(mc.cell_diameter_um, mc.c_internal, bins=[xb, yb], weights=z)
    N, _, _ = np.histogram2d(mc.cell_diameter_um, mc.c_internal, bins=[xb, yb])
    M = np.where(N > 0, H / np.maximum(N, 1), np.nan)
    pc = ax.pcolormesh(xb, yb, M.T, cmap="Blues", shading="flat",
                       norm=LogNorm(vmin=np.nanmin(M[M > 0]), vmax=np.nanmax(M)))
    cb = fig.colorbar(pc, ax=ax, pad=0.02)
    cb.set_label("mean " + OUTPUT_LABELS.get(surface_output, surface_output), fontsize=7)
    ax.set_yscale("log")
    ax.set_xlabel(r"Cell diameter $d$ ($\mu$m)")
    ax.set_ylabel(r"Pigment conc. $c_i$ (kg m$^{-3}$)")
    ax.set_title("MC interaction surface", loc="left")
    ax.grid(False)
    _tag(ax, "b")

    # (c) second-order indices
    ax = fig.add_subplot(gs[2])
    if s2 is not None:
        S = s2.copy().astype(float)
        S = S.where(~np.isnan(S), S.T)               # SALib fills the upper triangle only
        keep = [n for n in (s2_params or S2_DEFAULT) if n in S.index]
        S = S.loc[keep, keep].clip(lower=0)
        names = list(S.index)
        mask = np.tril(np.ones_like(S, dtype=bool))
        A = np.ma.array(S.to_numpy(), mask=mask)
        im = ax.imshow(A, cmap="Blues", vmin=0, vmax=max(0.05, float(np.nanmax(S.to_numpy()))))
        ax.set_xticks(range(len(names)))
        ax.set_yticks(range(len(names)))
        ax.set_xticklabels([SHORT.get(n, n) for n in names], rotation=45, ha="right")
        ax.set_yticklabels([SHORT.get(n, n) for n in names])
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                v = S.iloc[i, j]
                if v >= 0.01:
                    ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center", va="center", fontsize=6,
                            color="white" if v > 0.6 * np.nanmax(S.to_numpy()) else INK)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label(r"$S_{ij}$, " + OUTPUT_LABELS.get(s2_output, s2_output), fontsize=7)
        ax.grid(False)
        ax.set_title(r"Pairwise $S_{ij}$", loc="left")
    _tag(ax, "c")
    if demo:
        fig.suptitle("DEMO input - not Phase 1 results", color="#b00020")
    return fig


def fig_s3_convergence(conv: pd.DataFrame, space, output: str):
    meta = {p.name: p for p in space.params}
    fig, ax = plt.subplots(figsize=(F.SINGLE_COL + 0.9, 2.8), constrained_layout=True)
    top = conv[conv.n_base == conv.n_base.max()].sort_values("ST", ascending=False).parameter[:5]
    styles = ["-", "--", "-.", ":", (0, (5, 1, 1, 1))]
    for n, ls in zip(top, styles):
        c = conv[conv.parameter == n]
        col = GROUP_COLORS[meta[n].group] if n in meta else INK
        ax.errorbar(c.n_base, c.ST, yerr=c.ST_conf, color=col, ls=ls, marker="o", ms=3, capsize=2,
                    label=meta[n].label if n in meta else n)
    ax.set_xscale("log", base=2)
    ax.set_xlabel(r"Saltelli base sample size $N$")
    ax.set_ylabel(r"$S_T$ (95 % bootstrap CI)")
    ax.set_title(OUTPUT_LABELS.get(output, output), loc="left")
    ax.legend(fontsize=6.5)
    return fig
