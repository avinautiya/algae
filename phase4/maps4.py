"""
Figures 4A-4C (+ S4 validation, S5 diagnostics) for the Sentinel-2 inversion.

Maps are drawn in the scene's native UTM projection with cartopy (lat/lon graticule
labels, scale bar, north arrow); without cartopy they fall back to plain axes in
UTM kilometres. Colour use follows the earlier phases: one single-hue sequential map
per magnitude, a diverging map with a neutral grey midpoint for signed differences.
All figures: vector PDF + 600-dpi PNG via the Phase 2 style.
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import figures as F  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

INK, INK_2 = F.INK, F.INK_2
DIVERGING = LinearSegmentedColormap.from_list("bgr", ["#1f5fae", "#9ec5f0", "#e8e8e6", "#f2a07b", "#b33a0c"])
CMAPS = dict(biomass="Purples", sd="Oranges", albedo="Greys_r", rf="Reds", chi2="Greys")

try:
    import cartopy.crs as ccrs
    HAVE_CARTOPY = True
except ImportError:          # pragma: no cover
    HAVE_CARTOPY = False

set_style = F.set_style
save = F.save


class MapGrid:
    """Geometry helper: scene extent, projection and decorations."""

    def __init__(self, transform, crs, shape):
        self.t, self.crs, self.shape = transform, crs, shape
        H, W = shape
        self.extent = (transform.c, transform.c + W * transform.a, transform.f + H * transform.e, transform.f)
        self.proj = None
        if HAVE_CARTOPY:
            epsg = crs.to_epsg()
            if epsg and (32601 <= epsg <= 32660 or 32701 <= epsg <= 32760):
                self.proj = ccrs.UTM(epsg % 100, southern_hemisphere=epsg > 32700)

    def axes(self, fig, spec):
        if self.proj is not None:
            return fig.add_subplot(spec, projection=self.proj)
        return fig.add_subplot(spec)

    def show(self, ax, arr, cmap, norm=None, vmin=None, vmax=None, label="", scalebar=True, north=False,
             grid_labels=(True, True)):
        kw = dict(origin="upper", extent=self.extent, cmap=cmap, norm=norm, interpolation="nearest")
        if norm is None:
            kw.update(vmin=vmin, vmax=vmax)
        if self.proj is not None:
            im = ax.imshow(np.ma.masked_invalid(arr), transform=self.proj, **kw)
            ax.set_extent(self.extent, crs=self.proj)
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.3, color="#9a9a96",
                              alpha=0.7, x_inline=False, y_inline=False)
            gl.top_labels = gl.right_labels = False
            gl.bottom_labels, gl.left_labels = grid_labels
            gl.xlabel_style = gl.ylabel_style = dict(size=5.5, color=INK_2)
            gl.rotate_labels = False
        else:
            e = self.extent
            im = ax.imshow(np.ma.masked_invalid(arr), **dict(kw, extent=[v / 1000 for v in e]))
            ax.set_xlabel("UTM easting (km)")
            ax.set_ylabel("UTM northing (km)")
        ax.set_facecolor("#d9d9d6")
        cb = plt.colorbar(im, ax=ax, shrink=0.82, pad=0.03)
        cb.set_label(label, fontsize=6.5)
        cb.ax.tick_params(labelsize=6)
        if scalebar:
            self.scalebar(ax)
        if north:
            self.north_arrow(ax)
        return im

    def scalebar(self, ax, frac=0.25):
        x0, x1, y0, y1 = self.extent
        width = x1 - x0
        nice = [100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000]
        L = max(n for n in nice if n <= frac * width)
        s = 1.0 if self.proj is not None else 1e-3
        bx, by = x0 + 0.05 * width, y0 + 0.06 * (y1 - y0)
        h = 0.018 * (y1 - y0)
        kw = dict(transform=self.proj) if self.proj is not None else {}
        ax.add_patch(plt.Rectangle((bx * s, by * s), L * s / 2, h * s, facecolor=INK, edgecolor=INK, lw=0.5, **kw))
        ax.add_patch(plt.Rectangle(((bx + L / 2) * s, by * s), L * s / 2, h * s, facecolor="white", edgecolor=INK,
                                   lw=0.5, **kw))
        txt = f"{L / 1000:g} km" if L >= 1000 else f"{L} m"
        ax.text((bx + L / 2) * s, (by + 1.6 * h) * s, txt, ha="center", va="bottom", fontsize=6, color=INK,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7), **kw)

    def north_arrow(self, ax):
        ax.annotate("N", xy=(0.93, 0.93), xytext=(0.93, 0.80), xycoords="axes fraction",
                    ha="center", va="center", fontsize=7, fontweight="bold", color=INK,
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.0))


def _tag(ax, t, title=""):
    ax.set_title(f"$\\bf{{{t}}}$  {title}", loc="left", fontsize=7.5)


# --------------------------------------------------------------------------- #
# Figure 4A: biomass                                                            #
# --------------------------------------------------------------------------- #
def fig4a(g: MapGrid, panels: dict, sd_map, f_map, truth=None, rgb=None, demo_note="", f_info=None):
    """panels: label -> log10 B map (empirical, Tier A Bayesian, ours). sd_map: ours 95 % CI width."""
    fig = plt.figure(figsize=(F.DOUBLE_COL, 5.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    vals = np.concatenate([v[np.isfinite(v)].ravel() for v in panels.values()] +
                          ([truth[np.isfinite(truth)].ravel()] if truth is not None else []))
    vmin, vmax = np.percentile(vals, [1, 99])
    lab = r"$\log_{10}B$ (cells mL$^{-1}$)"
    ax = g.axes(fig, gs[0, 0])
    if truth is not None:
        g.show(ax, truth, CMAPS["biomass"], vmin=vmin, vmax=vmax, label=lab, north=True, grid_labels=(False, True))
        _tag(ax, "a", "synthetic truth")
    elif rgb is not None:
        ax.imshow(rgb, origin="upper", extent=g.extent, transform=g.proj) if g.proj is not None else \
            ax.imshow(rgb, origin="upper", extent=[v / 1000 for v in g.extent])
        g.scalebar(ax)
        g.north_arrow(ax)
        _tag(ax, "a", "Sentinel-2 true colour (L2A)")
    for (name, arr), spec, t in zip(panels.items(), [gs[0, 1], gs[0, 2], gs[1, 0]], "bcd"):
        ax = g.axes(fig, spec)
        g.show(ax, arr, CMAPS["biomass"], vmin=vmin, vmax=vmax, label=lab, grid_labels=(t == "d", t == "d"))
        _tag(ax, t, name)
    ax = g.axes(fig, gs[1, 1])
    g.show(ax, sd_map, CMAPS["sd"], label=r"95 % CI width of $\log_{10}B$ (dex)", grid_labels=(True, False))
    _tag(ax, "e", "ours: 95 % CI width")
    ax = g.axes(fig, gs[1, 2])
    g.show(ax, f_map, "Greens", vmin=0, vmax=1, label=r"posterior mean $f_{nordenskioeldii}$",
           grid_labels=(True, False))
    _tag(ax, "f", "community fraction")
    if f_info is not None:
        ax.text(0.5, 0.5, f"posterior/prior SD = {f_info:.2f}" + ("\n(not identifiable from 4 bands)" if f_info > 0.8
                                                                  else ""),
                transform=ax.transAxes, ha="center", va="center", fontsize=6.5, color=INK,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))
    if demo_note:
        fig.suptitle(demo_note, color="#b00020", fontsize=8)
    return fig


# --------------------------------------------------------------------------- #
# Figure 4B: albedo and forcing                                                 #
# --------------------------------------------------------------------------- #
def fig4b(g: MapGrid, bba, rf, drf, bba_sd, rf_sd, diag, diag_label, demo_note=""):
    fig = plt.figure(figsize=(F.DOUBLE_COL, 5.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)
    specs = [
        (gs[0, 0], bba, CMAPS["albedo"], None, r"broadband albedo $\alpha_{300-2500}$", "a", "albedo", (False, True)),
        (gs[0, 1], rf, CMAPS["rf"], None, r"$RF_{algae}$ (W m$^{-2}$)", "b", "algal radiative forcing", (False, False)),
        (gs[0, 2], drf, DIVERGING, _div_norm(drf), r"$RF_{ours}-RF_{A}$ (W m$^{-2}$)", "c",
         "forcing anomaly vs. Tier A", (False, False)),
        (gs[1, 0], bba_sd, CMAPS["sd"], None, r"posterior SD of $\alpha$", "d", "albedo uncertainty", (True, True)),
        (gs[1, 1], rf_sd, CMAPS["sd"], None, r"posterior SD of $RF$ (W m$^{-2}$)", "e", "forcing uncertainty",
         (True, False)),
        (gs[1, 2], diag, "Greys", None, diag_label, "f", "model evidence", (True, False)),
    ]
    for k, (spec, arr, cmap, norm, lab, t, title, gl) in enumerate(specs):
        ax = g.axes(fig, spec)
        g.show(ax, arr, cmap, norm=norm, label=lab, north=(k == 0), grid_labels=gl)
        _tag(ax, t, title)
    if demo_note:
        fig.suptitle(demo_note, color="#b00020", fontsize=8)
    return fig


def _div_norm(a):
    v = np.nanpercentile(np.abs(a), 98) if np.isfinite(a).any() else 1.0
    return TwoSlopeNorm(vcenter=0.0, vmin=-max(v, 1e-6), vmax=max(v, 1e-6))


# --------------------------------------------------------------------------- #
# Figure 4C: corner plot                                                        #
# --------------------------------------------------------------------------- #
LABELS = {"log_b": r"$\log_{10}B$", "f_n": r"$f_n$", "r_um": r"$r$ ($\mu$m)", "dust_ppb": r"$\log_{10}$dust",
          "k": r"$k$", "log_pig": r"$\log_{10}$pigment"}


def corner(samples, names, grid_marginals=None, truth=None, diag=None, title=""):
    """Corner plot: 2-D histograms (MCMC), 1-D histograms with the exact grid marginals
    overlaid (lines), truth markers if known, and per-parameter R-hat / ESS."""
    n = len(names)
    fig, axes = plt.subplots(n, n, figsize=(F.SINGLE_COL + 1.6, F.SINGLE_COL + 1.6), constrained_layout=True)
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if j > i:
                ax.axis("off")
                continue
            xi = samples[:, j]
            if i == j:
                ax.hist(xi, bins=40, density=True, color=F.SEQ_BLUE[1], alpha=0.55, lw=0)
                q = np.percentile(xi, [2.5, 50, 97.5])
                for v, ls in zip(q, (":", "-", ":")):
                    ax.axvline(v, color=INK, lw=0.7, ls=ls)
                if grid_marginals and names[j] in grid_marginals:
                    gx, gp = grid_marginals[names[j]]
                    dx = np.gradient(gx)
                    ax.plot(gx, gp / dx, color=F.TIERS["C"][0], lw=1.2)
                if truth is not None and names[j] in truth:
                    ax.axvline(truth[names[j]], color=F.TIERS["B"][0], lw=1.2)
                ax.set_yticks([])
                t = f"{q[1]:.3g}$^{{+{q[2] - q[1]:.2g}}}_{{-{q[1] - q[0]:.2g}}}$"
                if diag is not None:
                    t += f"\n$\\hat R$={diag['rhat'][j]:.3f}, ESS={diag['ess'][j]:.0f}"
                ax.set_title(t, fontsize=5.5)
            else:
                yi = samples[:, i]
                ax.hist2d(xi, yi, bins=35, cmap="Blues", cmin=1)
                if truth is not None and names[j] in truth and names[i] in truth:
                    ax.plot(truth[names[j]], truth[names[i]], marker="x", color=F.TIERS["B"][0], ms=5, mew=1.2)
            ax.grid(False)
            if i == n - 1:
                ax.set_xlabel(LABELS.get(names[j], names[j]), fontsize=7)
            else:
                ax.set_xticklabels([])
            if j == 0 and i > 0:
                ax.set_ylabel(LABELS.get(names[i], names[i]), fontsize=7)
            elif j > 0:
                ax.set_yticklabels([])
            ax.tick_params(labelsize=5.5)
    if title:
        fig.suptitle(title, fontsize=7.5)
    return fig


# --------------------------------------------------------------------------- #
# Supplementary                                                                 #
# --------------------------------------------------------------------------- #
def fig_s4_validation(truth, ests: dict, metrics_df, quantity="log_b", label=r"$\log_{10}B$", title=""):
    """Retrieved vs true for each method, with bias/RMSE (and coverage) annotated."""
    fig, axes = plt.subplots(1, len(ests), figsize=(F.DOUBLE_COL, 2.6), constrained_layout=True, sharey=True)
    allv = np.concatenate([np.ravel(truth)] + [np.ravel(np.asarray(e, float)) for e in ests.values()])
    lo, hi = np.nanpercentile(allv, [0.5, 99.5])
    pad = 0.05 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    cols = [F.TIERS["A"][0], F.SEQ_BLUE[2], F.TIERS["C"][0], F.TIERS["D"][0]]
    for ax, (name, e), col in zip(np.atleast_1d(axes), ests.items(), cols):
        if np.size(truth) > 200:
            ax.hexbin(np.ravel(truth), np.ravel(e), gridsize=45, cmap="Greys", mincnt=1, extent=(lo, hi, lo, hi))
        else:
            ax.scatter(np.ravel(truth), np.ravel(e), s=14, color=col, edgecolor="white", lw=0.4, zorder=3)
        ax.plot([lo, hi], [lo, hi], color=col, lw=1.2)
        m = metrics_df[(metrics_df.method == name) & ((metrics_df.quantity == quantity)
                                                      if "quantity" in metrics_df else True)]
        if len(m):
            m = m.iloc[0]
            t = f"bias {m.bias:+.2f}\nRMSE {m.rmse:.2f}\n$R^2$ {m.r2:.2f}"
            if "coverage95" in m and np.isfinite(m.get("coverage95", np.nan)):
                t += f"\n95% cov. {100 * m.coverage95:.0f}%"
            ax.text(0.04, 0.96, t, transform=ax.transAxes, va="top", fontsize=6.5)
        ax.set_title(name, loc="left", fontsize=7.5)
        ax.set_xlabel("true " + label)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    np.atleast_1d(axes)[0].set_ylabel("retrieved " + label)
    if title:
        fig.suptitle(title, fontsize=8)
    return fig
