"""
Publication figures for Phase 2 (Figures 2A-2C plus supplementary S1-S2).

Design rules applied throughout:
  * one y-scale per panel (no dual axes); related quantities go in adjacent panels;
  * each model tier keeps one fixed colour in every figure (validated
    colour-blind-safe categorical set), plus a distinct marker and line style and
    a tier letter in the legend, so identity never relies on colour alone; lines
    are directly labelled where the labels do not collide (Fig. 2C);
  * magnitude (cell size) uses a single-hue sequential ramp;
  * text is set in neutral ink; mathtext labels by default, true LaTeX with usetex=True;
  * every figure is saved as vector PDF and 600-dpi PNG.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

INK = "#1f1f1f"
INK_2 = "#52514e"
GRID = "#d9d9d6"
TIERS = {
    #       colour      marker  linestyle  label
    "A": ("#2a78d6", "o", "-",  "A  BioSNICAR empirical"),
    "B": ("#eb6834", "s", "--", "B  raw molecular MAC"),
    "C": ("#1baf7a", "^", "-",  "C  packaging-corrected"),
    "D": ("#eda100", "D", "-.", "D  Fe-phenolic + aggregated"),
}
SEQ_BLUE = ["#9ec5f0", "#4f94e0", "#1f5fae", "#123a6b"]     # light -> dark, one hue
SINGLE_COL, DOUBLE_COL = 3.5, 7.2                             # inches


def set_style(usetex: bool = False, font_size: float = 8.5):
    try:
        import seaborn as sns
        sns.set_theme(style="ticks", context="paper")
    except ImportError:
        pass
    plt.rcParams.update({
        "text.usetex": usetex,
        "font.family": "serif" if usetex else "sans-serif",
        "font.size": font_size, "axes.labelsize": font_size, "axes.titlesize": font_size,
        "legend.fontsize": font_size - 1, "xtick.labelsize": font_size - 0.5,
        "ytick.labelsize": font_size - 0.5,
        "axes.edgecolor": INK_2, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK_2, "ytick.color": INK_2, "axes.linewidth": 0.8,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
        "lines.linewidth": 1.6, "lines.markersize": 4.5,
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
        "pdf.fonttype": 42, "ps.fonttype": 42,           # editable text in Illustrator
        "mathtext.default": "regular",
    })


def save(fig, outdir: str, name: str, dpi: int = 600):
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for ext in ("pdf", "png"):
        p = os.path.join(outdir, f"{name}.{ext}")
        fig.savefig(p, dpi=dpi)
        paths.append(p)
    plt.close(fig)
    return paths


def _panel_tag(ax, tag):
    ax.text(-0.13, 1.02, tag, transform=ax.transAxes, fontweight="bold", va="bottom", ha="left")


def _end_label(ax, x, y, text, dy=0.0):
    ax.annotate(text, (x[-1], y[-1]), xytext=(4, dy), textcoords="offset points",
                va="center", ha="left", fontsize=plt.rcParams["legend.fontsize"], color=INK)


# --------------------------------------------------------------------------- #
# Figure 2A                                                                     #
# --------------------------------------------------------------------------- #
def fig2a(wl_nm, mac_raw: dict, mac_vivo: dict, q_star: dict, outdir: str,
          empirical_q=None, c_internal=None, demo=False):
    """(a) raw molecular MAC vs packaging-corrected MAC_vivo; (b) packaging factor Q*.

    mac_raw  : {label: MAC array} (e.g. Level 1 / Level 2), in m^2 kg^-1
    mac_vivo : {cell label: MAC_vivo array} for the Level 2 pigment
    q_star   : {cell label: Q* array}
    empirical_q : optional (wl_nm, Q) of BioSNICAR's measured glacier-algae packaging factor
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.9), constrained_layout=True)
    scale = 1e-3
    raw_styles = [("-", INK), (":", INK_2)]
    for (lab, m), (ls, col) in zip(mac_raw.items(), raw_styles):
        ax1.plot(wl_nm, m * scale, ls=ls, color=col, lw=1.6, label=f"{lab} (solution)")
    for (lab, m), col in zip(mac_vivo.items(), SEQ_BLUE[1:] if len(mac_vivo) <= 3 else SEQ_BLUE):
        ax1.plot(wl_nm, m * scale, color=col, lw=1.6, label=rf"in vivo, {lab}")
    ax1.set_xlim(wl_nm.min(), wl_nm.max())
    ax1.set_ylim(bottom=0)
    ax1.set_xlabel(r"Wavelength $\lambda$ (nm)")
    ax1.set_ylabel(r"MAC$(\lambda)$ ($10^{3}\ \mathrm{m^{2}\,kg^{-1}}$)")
    ax1.legend(loc="upper right")
    _panel_tag(ax1, "a")

    for (lab, q), col in zip(q_star.items(), SEQ_BLUE[1:] if len(q_star) <= 3 else SEQ_BLUE):
        ax2.plot(wl_nm, q, color=col, lw=1.6, label=lab)
    if empirical_q is not None:
        ax2.plot(empirical_q[0], empirical_q[1], color=INK_2, ls=(0, (1, 1.5)), lw=1.4,
                 label="BioSNICAR empirical (glacier algae)")
    ax2.axhline(1.0, color=INK_2, lw=0.6)
    ax2.set_xlim(wl_nm.min(), wl_nm.max())
    ax2.set_ylim(0, 1.05)
    ax2.set_xlabel(r"Wavelength $\lambda$ (nm)")
    ax2.set_ylabel(r"Packaging factor $Q^{*}=\mathrm{MAC_{vivo}}/\mathrm{MAC}$")
    ax2.legend(loc="lower right")
    if c_internal is not None:
        ax2.set_title(rf"intracellular pigment $c_i = {c_internal:g}$ kg m$^{{-3}}$", loc="left")
    _panel_tag(ax2, "b")
    if demo:
        fig.suptitle("DEMO input (BioSNICAR ppg.csv) - not Phase 1 results", color="#b00020")
    return fig


# --------------------------------------------------------------------------- #
# Figure 2B                                                                     #
# --------------------------------------------------------------------------- #
def fig2b(df: pd.DataFrame, ref: dict, outdir: str, demo=False):
    """Broadband (300-2500 nm) albedo vs algal concentration, tiers A-D, at the reference
    ice / illumination state; shaded band = range over the density sweep."""
    fig, ax = plt.subplots(figsize=(SINGLE_COL + 0.9, 2.9), constrained_layout=True)
    sel = df[(df.grain_um == ref["grain_um"]) & (df.sza == ref["sza"])]
    clean = sel[(sel.tier == "clean") & (sel.rho_top == ref["rho_top"])].bba.iloc[0]
    ax.axhline(clean, color=INK_2, lw=1.0, ls=(0, (4, 2)))
    ax.text(df.conc.max(), clean + 0.006, "clean ice", color=INK_2, va="bottom", ha="right")
    for t, (col, mk, ls, lab) in TIERS.items():
        s = sel[sel.tier == t]
        if s.empty:
            continue
        band = s.groupby("conc").bba.agg(["min", "max"])
        ax.fill_between(band.index, band["min"], band["max"], color=col, alpha=0.15, lw=0)
        r = s[s.rho_top == ref["rho_top"]].sort_values("conc")
        ax.plot(r.conc, r.bba, color=col, ls=ls, marker=mk, markevery=2, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel(r"Algal cell concentration (cells mL$^{-1}$)")
    ax.set_ylabel(r"Broadband albedo $\alpha_{300-2500}$")
    ax.legend(loc="lower left")
    ax.set_title(rf"${_r_symbol(ref)}={ref['grain_um'] / 1000:g}$ mm, $\rho={ref['rho_top']:g}$ kg m$^{{-3}}$, "
                 rf"SZA$={ref['sza']:g}^\circ$" + ("  [DEMO]" if demo else ""), loc="left")
    return fig


# --------------------------------------------------------------------------- #
# Figure 2C                                                                     #
# --------------------------------------------------------------------------- #
def fig2c(df: pd.DataFrame, ref: dict, concs, outdir: str, baseline: str = "A", demo=False):
    """Radiative-forcing anomaly RF_tier - RF_baseline vs grain size (small multiples per
    concentration); shaded band = range over the density sweep."""
    concs = list(concs)
    fig, axes = plt.subplots(1, len(concs), figsize=(DOUBLE_COL, 2.8), constrained_layout=True,
                             squeeze=False)
    sel = df[df.sza == ref["sza"]]
    base = sel[sel.tier == baseline][["grain_um", "rho_top", "conc", "rf"]].rename(columns={"rf": "rf0"})
    m = sel.merge(base, on=["grain_um", "rho_top", "conc"])
    m["drf"] = m.rf - m.rf0
    for ax, c, tag in zip(axes[0], concs, "abcdef"):
        mc = m[np.isclose(m.conc, c)]
        for t, (col, mk, ls, lab) in TIERS.items():
            if t == baseline:
                continue
            s = mc[mc.tier == t]
            if s.empty:
                continue
            band = s.groupby("grain_um").drf.agg(["min", "max"])
            x = band.index.to_numpy() / 1000.0
            ax.fill_between(x, band["min"], band["max"], color=col, alpha=0.15, lw=0)
            r = s[s.rho_top == ref["rho_top"]].sort_values("grain_um")
            ax.plot(r.grain_um / 1000.0, r.drf, color=col, ls=ls, marker=mk, label=lab)
            _end_label(ax, (r.grain_um / 1000.0).to_numpy(), r.drf.to_numpy(), t)
        ax.axhline(0, color=INK_2, lw=0.8)
        ax.set_xlabel(_r_label(ref))
        ax.set_ylabel(rf"$\Delta RF = RF_{{X}} - RF_{{{baseline}}}$ (W m$^{{-2}}$)")
        ax.set_title(_conc_title(c) + rf", SZA$={ref['sza']:g}^\circ$", loc="left")
        _panel_tag(ax, tag)
    axes[0][0].legend(loc="best")
    if demo:
        fig.suptitle("DEMO input - not Phase 1 results", color="#b00020")
    return fig


def _r_symbol(ref):
    return "r_{bubble}" if ref.get("ice_mode") == "bubbly" else "r_{grain}"


def _r_label(ref):
    if ref.get("ice_mode") == "bubbly":
        return r"Bubbly-ice optical radius $r_{bubble}$ (mm)"
    return r"Ice grain radius $r_{grain}$ (mm)"


def _conc_title(c):
    e = int(np.floor(np.log10(c)))
    mant = c / 10**e
    mtxt = "" if np.isclose(mant, 1) else rf"{mant:g}\times"
    return rf"${mtxt}10^{{{e}}}$ cells mL$^{{-1}}$"


# --------------------------------------------------------------------------- #
# Supplementary                                                                 #
# --------------------------------------------------------------------------- #
def fig_s1_cell_optics(wl_nm, optics: dict, outdir: str):
    """Per-cell absorption cross-section and single-scattering albedo for tiers A-D."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.8), constrained_layout=True)
    for t, (col, mk, ls, lab) in TIERS.items():
        if t not in optics:
            continue
        o = optics[t]
        ax1.plot(wl_nm, o["abs_xsc"] * 1e12, color=col, ls=ls, label=lab)
        ax2.plot(wl_nm, o["ss_alb"], color=col, ls=ls, label=lab)
    for ax in (ax1, ax2):
        ax.set_xlim(300, 1000)
        ax.set_xlabel(r"Wavelength $\lambda$ (nm)")
    ax1.set_yscale("log")
    ax1.set_ylim(1e-1, None)
    ax1.set_ylabel(r"$\sigma_{abs}$ per cell ($\mu\mathrm{m^{2}}$)")
    ax2.set_ylabel(r"Single-scattering albedo $\omega$")
    ax2.set_ylim(0, 1.02)
    ax1.legend(loc="lower left")
    _panel_tag(ax1, "a")
    _panel_tag(ax2, "b")
    return fig


def fig_s2_spectral_albedo(wl_nm, albedo: dict, clean, conc, outdir: str):
    fig, ax = plt.subplots(figsize=(SINGLE_COL + 0.9, 2.8), constrained_layout=True)
    ax.plot(wl_nm, clean, color=INK_2, ls=(0, (4, 2)), lw=1.2, label="clean ice")
    for t, (col, mk, ls, lab) in TIERS.items():
        if t in albedo:
            ax.plot(wl_nm, albedo[t], color=col, ls=ls, label=lab)
    ax.set_xlim(300, 2500)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"Wavelength $\lambda$ (nm)")
    ax.set_ylabel(r"Spectral albedo $\alpha(\lambda)$")
    ax.set_title(_conc_title(conc), loc="left")
    ax.legend(loc="upper right")
    return fig
