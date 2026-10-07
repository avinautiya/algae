"""
Fast packaging factor Q* for Monte Carlo / Sobol work.

For a convex cell, Q* depends on only two dimensionless numbers:

    Q*(tau, shape) = < 1 - exp(-tau * l_hat) > / tau,   tau = a_i <l>,  l_hat = l / <l>

where <l> = 4V/S is the mean chord. For a cylinder the chord distribution of
l_hat depends only on the aspect ratio L/d. We therefore tabulate Q* once on a
(log aspect ratio x log tau) grid from the exact Monte Carlo chords of Phase 2
(pigment_packaging.chords_cylinder) and interpolate bilinearly. This replaces a
~1 s chord average per call with microseconds while matching the direct method
to <0.3 % (tests/test_phase3.py).
"""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy.interpolate import RegularGridInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

from pigment_packaging import CellGeometry, chords_cylinder, q_star_sphere  # noqa: E402


class QStarTable:
    def __init__(self, aspect_min: float = 0.5, aspect_max: float = 8.0, n_aspect: int = 41,
                 tau_min: float = 1e-4, tau_max: float = 1e5, n_tau: int = 361,
                 n_chords: int = 60_000, seed: int = 7, cache: str | None = None):
        self.log_aspect = np.linspace(np.log(aspect_min), np.log(aspect_max), n_aspect)
        self.log_tau = np.linspace(np.log(tau_min), np.log(tau_max), n_tau)
        if cache and os.path.isfile(cache):
            d = np.load(cache)
            if d["table"].shape == (n_aspect, n_tau) and np.allclose(d["log_aspect"], self.log_aspect):
                self.table = d["table"]
                self._build_interp()
                return
        tau = np.exp(self.log_tau)
        table = np.empty((n_aspect, n_tau))
        for i, la in enumerate(self.log_aspect):
            ar = np.exp(la)                               # L/d
            r = 1.0
            l = chords_cylinder(r, 2.0 * r * ar, n=n_chords, seed=seed + i)
            lh = l / l.mean()
            # chunk over tau to bound memory
            for j0 in range(0, n_tau, 40):
                tt = tau[j0:j0 + 40, None]
                table[i, j0:j0 + 40] = (-np.expm1(-tt * lh[None, :])).mean(axis=1) / tt[:, 0]
        self.table = table
        if cache:
            np.savez(cache, table=table, log_aspect=self.log_aspect, log_tau=self.log_tau)
        self._build_interp()

    def _build_interp(self):
        self._interp = RegularGridInterpolator((self.log_aspect, self.log_tau), np.log(self.table),
                                               bounds_error=False, fill_value=None)

    def __call__(self, a_internal_per_m, geom: CellGeometry):
        """Drop-in replacement for pigment_packaging.q_star(a, geom)."""
        a = np.asarray(a_internal_per_m, dtype=float)
        if geom.shape == "sphere":
            return q_star_sphere(a * 2.0 * geom.radius * 1e-6)
        tau = a * geom.mean_chord_um * 1e-6
        out = np.ones_like(tau)
        pos = tau > np.exp(self.log_tau[0])
        la = np.clip(np.log(geom.length / (2.0 * geom.radius)), self.log_aspect[0], self.log_aspect[-1])
        lt = np.log(tau[pos])
        big = lt > self.log_tau[-1]
        q = np.empty(lt.size)
        pts = np.column_stack([np.full((~big).sum(), la), lt[~big]])
        q[~big] = np.exp(self._interp(pts))
        q[big] = 1.0 / tau[pos][big]                      # optically black cell: Q* -> 1/tau
        out[pos] = q
        return out
