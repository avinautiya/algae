"""Focused checks of units, spectral operators and parameter bounds (inversion_m1/operators.py)."""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import operators as O  # noqa: E402


# ------------------------------------------------------------------------------------------------ units
def test_meltwater_to_biosnicar_matches_the_bridge_constant():
    assert O.cells_meltwater_to_biosnicar(1000.0) == pytest.approx(917.0)
    assert O.cells_biosnicar_to_meltwater(O.cells_meltwater_to_biosnicar(3.7e4)) == pytest.approx(3.7e4)
    # the same constant as the earlier phases: a silent change would shift every abundance by 8 %
    src = open(os.path.join(HERE, "..", "..", "phase2", "biosnicar_bridge.py")).read()
    assert "MELTWATER_TO_BIOSNICAR = 917.0 / 1000.0" in src


def test_dust_ppb_ug_per_g_round_trip_and_value():
    assert O.ug_per_g_to_ppb(342.0) == pytest.approx(3.42e5)            # the S6 mean dust loading
    assert O.ppb_to_ug_per_g(O.ug_per_g_to_ppb(5.0)) == pytest.approx(5.0)


# ------------------------------------------------------------------------------ band-average operator
def test_band_average_constant_and_linear():
    wl = np.arange(350.0, 1501.0, 1.0)
    c = np.arange(400.0, 1300.1, 10.0)
    assert np.allclose(O.band_average(wl, np.full(wl.size, 0.37), c, 10.0), 0.37)
    # linear spectrum: the box average over [c-5, c+5) of integer-nm samples is c - 0.5
    y = 0.001 * (wl - 350.0)
    got = O.band_average(wl, y, c, 10.0)
    assert np.allclose(got, 0.001 * (c - 0.5 - 350.0))


def test_band_average_ignores_nan_and_never_extrapolates():
    wl = np.arange(350.0, 1000.0, 1.0)
    y = np.ones(wl.size)
    y[(wl >= 600) & (wl < 604)] = np.nan                                  # 6 of 10 finite: kept (>= 50 %)
    y[(wl >= 700) & (wl < 708)] = np.nan                                  # 2 of 10 finite: dropped
    got = O.band_average(wl, y, np.array([605.0, 705.0, 345.0, 1200.0]), 10.0)
    assert got[0] == pytest.approx(1.0) and np.isnan(got[1])
    assert np.isnan(got[2]) and np.isnan(got[3])                           # outside the observed range


def test_band_average_rejects_unsorted_wavelengths():
    with pytest.raises(ValueError):
        O.band_average(np.array([2.0, 1.0, 3.0]), np.ones(3), [2.0], 1.0)


# --------------------------------------------------------------------------------- SRF band operator
def test_srf_average_boxcar_and_scale_invariance():
    wl = np.arange(400.0, 900.0, 1.0)
    y = (wl - 400.0) / 500.0
    sw = np.arange(600.0, 701.0, 1.0)
    s = np.ones(sw.size)
    v = O.srf_average(wl, y, sw, s)
    assert v == pytest.approx((650.0 - 400.0) / 500.0, abs=1e-9)           # mean of a linear y over 600-700 nm = value at 650
    assert O.srf_average(wl, y, sw, 7.3 * s) == pytest.approx(v)           # normalised: SRF scale cancels


def test_srf_average_requires_coverage_and_nonnegative_response():
    wl = np.arange(500.0, 650.0, 1.0)
    sw = np.arange(600.0, 701.0, 1.0)
    with pytest.raises(ValueError):
        O.srf_average(wl, np.ones(wl.size), sw, np.ones(sw.size))          # response runs beyond the spectrum
    with pytest.raises(ValueError):
        O.srf_average(np.arange(500.0, 800.0), np.ones(300), sw, -np.ones(sw.size))


def test_srf_average_on_the_real_sentinel2_response():
    import pandas as pd
    srf = pd.read_csv(os.path.join(HERE, "..", "..", "data", "empirical", "esa_s2_srf_TN-15-0007_v4.0.csv"))
    wl = np.arange(350.0, 2500.0, 1.0)
    for band in ("S2A_B2", "S2A_B3", "S2A_B4", "S2A_B8"):
        v = O.srf_average(wl, np.full(wl.size, 0.42), srf.wavelength_nm.to_numpy(), srf[band].to_numpy())
        assert v == pytest.approx(0.42)                                    # a constant albedo is band-invariant


# ------------------------------------------------------------------------------------------ broadband
def test_broadband_constant_and_weighting():
    wl = np.arange(350.0, 2500.0, 1.0)
    fw = np.arange(300.0, 2501.0, 1.0)
    flux = np.exp(-((fw - 700.0) / 300.0) ** 2)
    assert O.broadband(wl, np.full(wl.size, 0.3), fw, flux) == pytest.approx(0.3)
    assert O.broadband(wl, np.full(wl.size, 0.3), fw, 5.0 * flux) == pytest.approx(0.3)      # flux scale cancels
    step = np.where(wl < 1000.0, 0.8, 0.1)                                                   # bright VIS, dark NIR
    b = O.broadband(wl, step, fw, flux)
    g = np.arange(300.0, 2501.0, 1.0)
    h = np.interp(g, wl, step)
    assert b == pytest.approx(np.sum(h * flux) / np.sum(flux))
    assert 0.1 < b < 0.8


def test_broadband_holds_the_end_value_outside_the_observed_range():
    wl = np.arange(350.0, 2500.0, 1.0)
    fw = np.arange(300.0, 2501.0, 1.0)
    flux = np.ones(fw.size)
    sp = np.where(wl <= 350.0, 0.2, 0.6)                                                      # 350 nm value is 0.2
    got = O.broadband(wl, sp, fw, flux)
    ref = (51 * 0.2 + 0.6 * (2201 - 51)) / 2201                                               # 300-349 held at 0.2
    assert got == pytest.approx(ref, abs=2e-3)


def test_broadband_rejects_bad_flux_and_handles_all_nan():
    wl = np.arange(350.0, 2500.0, 1.0)
    with pytest.raises(ValueError):
        O.broadband(wl, np.full(wl.size, 0.3), np.arange(300.0, 2501.0), -np.ones(2201))
    assert np.isnan(O.broadband(wl, np.full(wl.size, np.nan), np.arange(300.0, 2501.0), np.ones(2201)))


# --------------------------------------------------------------------- HCRF is not albedo, and bounds
def test_hcrf_albedo_round_trip_and_positive_k():
    a = np.array([0.2, 0.5, 0.8])
    assert np.allclose(O.albedo_from_hcrf(O.hcrf_from_albedo(a, 0.9), 0.9), a)
    with pytest.raises(ValueError):
        O.albedo_from_hcrf(a, 0.0)
    # HCRF cannot be passed off as albedo: the same number maps to different albedos for different k
    assert not np.allclose(O.albedo_from_hcrf(a, 0.8), O.albedo_from_hcrf(a, 1.1))


def test_check_albedo_bounds_flags_unmasked_absorption_bands():
    ok = np.array([0.0, 0.4, 0.95, np.nan])
    assert np.allclose(O.check_albedo_bounds(ok, name="x")[:3], ok[:3])
    for bad in (np.array([-308.0, 0.3]), np.array([0.3, 166.0]), np.array([0.3, 1.2])):
        with pytest.raises(ValueError):
            O.check_albedo_bounds(bad)


def test_window_mask_includes_and_excludes():
    wl = np.arange(300.0, 1500.0, 10.0)
    m = O.window_mask(wl, 400.0, 1300.0, exclude=[(750.0, 770.0)])
    assert wl[m].min() == 400.0 and wl[m].max() == 1300.0
    assert not m[(wl >= 750.0) & (wl <= 770.0)].any()


def test_params_from_config_validate_intervals():
    ok = O.make_params({"a": {"bounds": [0.0, 2.0], "unit": "x"}, "b": {"bounds": [-1.0, 3.0]}})
    assert [p.name for p in ok] == ["a", "b"]
    for bad in ({"a": {"bounds": [2.0, 1.0]}}, {"a": {"bounds": [0.0, float("inf")]}}, {"a": {"bounds": [1.0, 1.0]}}):
        with pytest.raises(ValueError):
            O.make_params(bad)


def test_param_unit_round_trip_clip_and_at_bound():
    p = O.Param("logB", 1.0, 6.5, "log10 cells/mL")
    assert p.from_unit(p.to_unit(3.2)) == pytest.approx(3.2)
    assert p.contains(1.0) and not p.contains(6.6)
    assert O.clip_to_bounds([0.0, 9.0], [p, p]).tolist() == [1.0, 6.5]
    assert p.at_bound(1.02) and p.at_bound(6.49) and not p.at_bound(3.0)
