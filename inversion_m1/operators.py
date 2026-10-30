"""
Units, spectral operators and parameter bounds for the glacier-darkening inversion study.

No radiative-transfer or file-format dependencies: everything here is plain NumPy so that the conventions
can be tested on their own (inversion_m1/tests/test_operators.py).

Conventions
-----------
* Wavelengths are in nm everywhere in this module (BioSNICAR's own grid is in um; convert at the boundary).
* Algal abundance is cells per mL of MELTWATER (how counts are reported). BioSNICAR converts its cells/mL input
  to cells per kg of ice as conc / 917 * 1e6, i.e. per mL of SOLID ice, so a meltwater count c is passed as
  0.917 * c (same convention as phase2/biosnicar_bridge.MELTWATER_TO_BIOSNICAR).
* Mineral dust and black carbon are mass mixing ratios in ppb (ng of particles per g of ice).
  1 ug/g = 1000 ppb.
* Spectral albedo and HCRF are dimensionless. HCRF is a DIRECTIONAL quantity; albedo is hemispherical. They are
  related only through an anisotropy factor k = HCRF / albedo, which is itself measured/uncertain.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MELTWATER_TO_BIOSNICAR = 917.0 / 1000.0


# ----------------------------------------------------------------------------------------------- units
def cells_meltwater_to_biosnicar(c):
    """Cells per mL of meltwater -> BioSNICAR's cells-per-mL input."""
    return np.asarray(c, dtype=float) * MELTWATER_TO_BIOSNICAR


def cells_biosnicar_to_meltwater(c):
    return np.asarray(c, dtype=float) / MELTWATER_TO_BIOSNICAR


def ug_per_g_to_ppb(x):
    return np.asarray(x, dtype=float) * 1.0e3


def ppb_to_ug_per_g(x):
    return np.asarray(x, dtype=float) / 1.0e3


# ----------------------------------------------------------------------------------- spectral operators
def band_average(wl_nm, y, centres_nm, width_nm: float, min_finite_fraction: float = 0.5):
    """Box average of a finely sampled spectrum onto coarser bands [c - w/2, c + w/2).

    This is the observation operator that maps a 1-nm field spectrum onto the model's 10-nm bands. Non-finite
    samples are ignored; a band with fewer than `min_finite_fraction` of its expected samples finite, or lying
    outside the observed range, is returned as NaN (never extrapolated).
    """
    wl = np.asarray(wl_nm, dtype=float)
    y = np.asarray(y, dtype=float)
    if wl.ndim != 1 or y.shape != wl.shape:
        raise ValueError("wl_nm and y must be 1-D arrays of the same length")
    if np.any(np.diff(wl) <= 0):
        raise ValueError("wavelengths must be strictly increasing")
    step = float(np.median(np.diff(wl)))
    expected = max(1.0, width_nm / step)
    out = np.full(len(centres_nm), np.nan)
    for j, c in enumerate(np.asarray(centres_nm, dtype=float)):
        lo, hi = c - width_nm / 2.0, c + width_nm / 2.0
        if lo < wl[0] - 1e-9 or hi > wl[-1] + step + 1e-9:
            continue
        m = (wl >= lo) & (wl < hi)
        v = y[m]
        v = v[np.isfinite(v)]
        if v.size >= max(1.0, min_finite_fraction * expected):
            out[j] = v.mean()
    return out


def srf_average(wl_nm, y, srf_wl_nm, srf):
    """Spectral-response-weighted band value: integral(y s) / integral(s), trapezoid on the SRF support.

    The spectrum must cover the whole support where the response is non-zero; otherwise the band value would be
    biased, so this raises instead of renormalising over a truncated response.
    """
    wl = np.asarray(wl_nm, dtype=float)
    y = np.asarray(y, dtype=float)
    sw = np.asarray(srf_wl_nm, dtype=float)
    s = np.asarray(srf, dtype=float)
    if np.any(s < 0):
        raise ValueError("spectral response must be non-negative")
    nz = s > 0
    if not nz.any():
        raise ValueError("spectral response is identically zero")
    lo, hi = sw[nz].min(), sw[nz].max()
    if lo < wl.min() - 1e-9 or hi > wl.max() + 1e-9:
        raise ValueError(f"spectrum ({wl.min():.0f}-{wl.max():.0f} nm) does not cover the response ({lo:.0f}-{hi:.0f} nm)")
    ok = np.isfinite(y)
    yi = np.interp(sw, wl[ok], y[ok])
    return float(np.trapezoid(yi * s, sw) / np.trapezoid(s, sw))


def broadband(wl_nm, spectrum, flux_wl_nm, flux, lo_nm: float = 300.0, hi_nm: float = 2500.0, step_nm: float = 1.0):
    """Irradiance-weighted broadband albedo, with the convention of phase4/h1_albedo.broadband: the finite
    spectral values are interpolated onto a 1-nm grid over [lo, hi] (end values held outside the observed range,
    e.g. 300-350 nm held at the 350 nm value) and weighted by the model's incoming spectral flux."""
    wl = np.asarray(wl_nm, dtype=float)
    sp = np.asarray(spectrum, dtype=float)
    ok = np.isfinite(sp)
    if ok.sum() < 2:
        return float("nan")
    grid = np.arange(lo_nm, hi_nm + 0.5 * step_nm, step_nm)
    h = np.interp(grid, wl[ok], sp[ok])
    f = np.interp(grid, np.asarray(flux_wl_nm, dtype=float), np.asarray(flux, dtype=float))
    if not np.all(f >= 0) or f.sum() <= 0:
        raise ValueError("incoming flux must be non-negative with a positive total")
    return float(np.sum(h * f) / np.sum(f))


def window_mask(wl_nm, lo_nm: float, hi_nm: float, exclude=()):
    """Boolean mask of wavelengths in [lo, hi] and outside every excluded (lo, hi) window."""
    wl = np.asarray(wl_nm, dtype=float)
    m = (wl >= lo_nm) & (wl <= hi_nm)
    for a, b in exclude:
        m &= ~((wl >= a) & (wl <= b))
    return m


def hcrf_from_albedo(albedo, k):
    """Directional reflectance factor implied by an albedo and an anisotropy factor k = HCRF / albedo."""
    return np.asarray(k, dtype=float) * np.asarray(albedo, dtype=float)


def albedo_from_hcrf(hcrf, k):
    k = np.asarray(k, dtype=float)
    if np.any(k <= 0):
        raise ValueError("anisotropy factor k must be positive")
    return np.asarray(hcrf, dtype=float) / k


def check_albedo_bounds(a, tol: float = 0.02, name: str = "albedo"):
    """Raise if any finite value lies outside [-tol, 1 + tol]: such values are non-physical for an albedo and
    signal an unmasked absorption band, a unit error or an HCRF passed as albedo."""
    a = np.asarray(a, dtype=float)
    f = a[np.isfinite(a)]
    if f.size and (f.min() < -tol or f.max() > 1.0 + tol):
        raise ValueError(f"{name} outside [0, 1] beyond tolerance: min {f.min():.3f}, max {f.max():.3f}")
    return a


# ----------------------------------------------------------------------------------- parameter bounds
@dataclass(frozen=True)
class Param:
    name: str
    lo: float
    hi: float
    unit: str

    def contains(self, x: float) -> bool:
        return self.lo <= x <= self.hi

    def to_unit(self, x):
        return (np.asarray(x, dtype=float) - self.lo) / (self.hi - self.lo)

    def from_unit(self, u):
        return self.lo + np.asarray(u, dtype=float) * (self.hi - self.lo)

    def at_bound(self, x: float, frac: float = 0.01) -> bool:
        w = self.hi - self.lo
        return bool(x <= self.lo + frac * w or x >= self.hi - frac * w)


def make_params(spec: dict) -> list[Param]:
    """Parameter list from the config's `parameters` block, validating that every interval is proper."""
    out = []
    for name, d in spec.items():
        lo, hi = (float(v) for v in d["bounds"])
        if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
            raise ValueError(f"parameter {name}: bounds must be finite and increasing, got {d['bounds']}")
        out.append(Param(name, lo, hi, d.get("unit", "")))
    return out


def clip_to_bounds(theta, params: list[Param]):
    return np.array([min(max(t, p.lo), p.hi) for t, p in zip(theta, params)], dtype=float)
