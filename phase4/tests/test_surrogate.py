"""Surrogate API: interpolation, SZA blending, domain flags."""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))

import emulator as E  # noqa: E402
import surrogate_api as SA  # noqa: E402


def _em(sza):
    axes = {"log_b": np.arange(1.0, 6.01, 0.5), "f_n": np.linspace(0, 1, 3), "r_um": np.array([1000.0, 4000.0]),
            "dust_ppb": np.array([3e4, 3e5])}
    g = np.meshgrid(*axes.values(), indexing="ij")
    bba = 0.65 - 0.03 * g[0] - 0.01 * g[1] - 1e-6 * g[2] - 0.001 * sza
    return E.Emulator(axes, {"bba": bba, "rf_algae": 10 * g[0] + sza, "bands": np.stack([bba] * 4, -1)}, {})


def test_interpolation_sza_blend_and_flags():
    s = SA.Surrogate({40.0: _em(40), 60.0: _em(60)}, "test")
    r = s(3.25, 50.0, r_um=2500.0, f_n=0.5, dust_ppb=1e5)
    assert r["in_domain"] and not r["flags"]
    assert np.isclose(r["bba"], 0.65 - 0.03 * 3.25 - 0.005 - 1e-6 * 2500 - 0.05)
    assert np.isclose(r["rf_algae"], 32.5 + 50.0)
    assert np.isclose(r["dalpha_algae_vs_lowest_abundance"], 0.03 * (3.25 - 1.0))
    out = s(7.5, 70.0)
    assert not out["in_domain"] and len(out["flags"]) == 3        # log_b, sza and the default dust 0 ppb
    assert np.isclose(out["rf_algae"], 60 + 60)                    # clipped to the boundary, flagged
