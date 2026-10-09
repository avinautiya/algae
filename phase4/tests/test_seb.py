"""Surface energy balance: longwave accounting, closure, time axis and gaps, missing shortwave, albedo
bounds, timestep convergence, paired runs."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))

import seb  # noqa: E402


def test_melting_surface_energy_closure_with_reflected_longwave():
    sw, lw, Ta, qa, U, p = 600.0, 300.0, 278.0, 0.004, 5.0, 870.0
    melt, Ts, H, LE = seb.solve_hour([sw], [lw], [Ta], [qa], [U], [p], [2.5])
    assert Ts[0] == seb.T0 and melt[0] > 0
    # absorbed LW is eps LW_down (the rest is reflected), emitted eps sigma T^4
    lw_net = seb.EPS_ICE * (lw - seb.SIGMA_SB * seb.T0 ** 4)
    assert np.isclose(melt[0], sw + lw_net + H[0] + LE[0])
    # the former accounting (all LW_down absorbed) differs by (1 - eps) LW_down = 6 W m-2 here
    assert abs((sw + lw - seb.EPS_ICE * seb.SIGMA_SB * seb.T0 ** 4 + H[0] + LE[0]) - melt[0] - 0.02 * lw) < 1e-9


def test_cold_surface_zero_net_flux_and_no_melt():
    melt, Ts, H, LE = seb.solve_hour([0.0], [220.0], [263.0], [0.0015], [3.0], [870.0], [2.5])
    assert melt[0] == 0 and Ts[0] < seb.T0
    q, _, _ = seb.balance(Ts, 0.0, 220.0, 263.0, 0.0015, 3.0, 870.0, 2.5, 1e-3)
    assert abs(q[0]) < 1e-6


def _diurnal(n=96):
    t = np.arange(n)
    sw = np.clip(700 * np.sin(2 * np.pi * (t - 6) / 24), 0, None)
    return (sw * 0.5, np.full(n, 270.0), 272.0 + 3 * np.sin(2 * np.pi * (t - 9) / 24), np.full(n, 0.0035),
            np.full(n, 4.0), np.full(n, 870.0), np.full(n, 2.5))


def test_energy_closure_and_timestep_convergence():
    args = _diurnal()
    totals = {}
    for n_sub in (1, 4, 16, 64):
        m, Ts, H, LE, b = seb.solve_series(*args, slab_m=0.1, n_sub=n_sub, return_budget=True)
        assert abs(b["closure_residual"]) < 1e-6 * max(1.0, abs(b["energy_in"]))   # exact bookkeeping
        totals[n_sub] = m.sum()
    rel = {k: abs(v - totals[64]) / totals[64] for k, v in totals.items()}
    assert rel[16] < 0.002 and rel[4] < 0.01                                        # converged sub-stepping


def test_slab_limits():
    args = _diurnal(48)
    m0, *_ = seb.solve_series(*args, slab_m=0.0)
    mi, *_ = seb.solve_hour(*args)
    assert np.allclose(m0, mi)
    m1, Ts1, *_ = seb.solve_series(*args, slab_m=0.1)
    assert m1.sum() < m0.sum()


def _forcing(times, sw=500.0):
    n = len(times)
    return pd.DataFrame(dict(time=pd.to_datetime(times), sw_down=np.full(n, sw), dlr=270.0, T_a=274.0, q_a=0.004,
                             wspd_u=4.0, p_u=870.0, z_meas=2.5))


def test_time_axis_gaps_and_missing_shortwave_are_kept():
    times = list(pd.date_range("2019-07-01", periods=6, freq="h")) + list(pd.date_range("2019-07-01 12:00", periods=4, freq="h"))
    f = _forcing(times)
    f.loc[2, "sw_down"] = np.nan                       # missing SW must NOT become zero
    r = seb.run(f, 0.5)
    assert len(r) == 16 and r.time.diff().dropna().eq(pd.Timedelta("1h")).all()     # actual hourly axis
    assert np.isnan(r.melt_w[2]) and not r.valid[2]
    assert (~r.valid).sum() == 1 + 6 and r.restart.sum() == 2                        # gap rows + restarts
    f2 = f.copy()
    f2.loc[3, "dlr"] = np.nan                          # short non-SW gap: interpolated, recorded
    assert seb.prepare(f2)[0].attrs["filled_hours"]["dlr"] == 1


def test_albedo_bounds_enforced():
    f = _forcing(pd.date_range("2019-07-01", periods=3, freq="h"))
    for bad in (1.2, -0.1, np.nan):
        with pytest.raises(ValueError):
            seb.run(f, bad)


def test_paired_runs_zero_without_difference_and_bounded_by_potential():
    f = _forcing(pd.date_range("2019-07-01", periods=72, freq="h"))
    f["sw_down"] = np.clip(700 * np.sin(2 * np.pi * (np.arange(72) - 6) / 24), 0, None)
    r0 = seb.paired_algae(2019, 0.5, 0.5, forcing=f)
    assert r0["modelled_algal_melt_increment_mwe"] == 0.0
    r1 = seb.paired_algae(2019, 0.47, 0.5, forcing=f)
    assert 0 < r1["modelled_algal_melt_increment_mwe"] <= 1.0001 * r1["potential_algal_melt_mwe"]
