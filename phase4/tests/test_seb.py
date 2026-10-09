"""Surface energy balance: closure, cold-surface solution, slab limit, paired runs."""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))

import seb  # noqa: E402


def test_melting_surface_energy_closure():
    sw, lw, Ta, qa, U, p = 600.0, 300.0, 278.0, 0.004, 5.0, 870.0
    melt, Ts, H, LE = seb.solve_hour([sw], [lw], [Ta], [qa], [U], [p], [2.5])
    assert Ts[0] == seb.T0 and melt[0] > 0
    lw_out = seb.EPS_ICE * seb.SIGMA_SB * seb.T0 ** 4
    assert np.isclose(melt[0], sw + lw - lw_out + H[0] + LE[0])


def test_cold_surface_zero_net_flux_and_no_melt():
    melt, Ts, H, LE = seb.solve_hour([0.0], [220.0], [263.0], [0.0015], [3.0], [870.0], [2.5])
    assert melt[0] == 0 and Ts[0] < seb.T0
    q, _, _ = seb.balance(Ts, 0.0, 220.0, 263.0, 0.0015, 3.0, 870.0, 2.5, 1e-3)
    assert abs(q[0]) < 1e-6


def test_slab_carries_cold_content_and_zero_slab_matches_instantaneous():
    n = 48
    t = np.arange(n)
    sw = np.clip(700 * np.sin(2 * np.pi * (t - 6) / 24), 0, None)
    args = (sw * 0.5, np.full(n, 270.0), 272.0 + 3 * np.sin(2 * np.pi * (t - 9) / 24), np.full(n, 0.0035),
            np.full(n, 4.0), np.full(n, 870.0), np.full(n, 2.5))
    m0, *_ = seb.solve_series(*args, slab_m=0.0)
    mi, *_ = seb.solve_hour(*args)
    assert np.allclose(m0, mi)
    m1, Ts1, *_ = seb.solve_series(*args, slab_m=0.1)
    assert m1.sum() < m0.sum()                              # night-time cold content delays morning melt
    assert Ts1.min() > seb.solve_hour(*args)[1].min()       # slab damps the night-time cooling


def test_paired_runs_zero_without_algae_and_positive_with():
    r0 = seb.paired_algae(2017, 0.55, 0.0, slab_m=0.1)
    assert r0["algal_melt_mwe"] == 0.0 and r0["potential_algal_melt_mwe"] == 0.0
    r1 = seb.paired_algae(2017, 0.55, 0.03, slab_m=0.1)
    assert 0 < r1["algal_melt_mwe"] <= 1.05 * r1["potential_algal_melt_mwe"]
