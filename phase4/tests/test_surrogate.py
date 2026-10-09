"""Surrogate API: interpolation, SZA blending, zero-algae reference, domain flags, qualification."""

import os
import sys

import numpy as np
import pytest

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
    sw = 800.0 - sza
    return E.Emulator(axes, {"bba": bba, "rf_algae": sw * 0.03 * g[0], "bands": np.stack([bba] * 4, -1)},
                      {"tag": f"t{sza}", "sw_down": sw})


def _ems():
    return {40.0: _em(40), 60.0: _em(60)}


def test_qualification_required():
    with pytest.raises(SA.QualificationError):
        SA.Surrogate(_ems(), "test")
    bad = dict(passed=True, emulator_tags={"40.0": "other", "60.0": "t60"})
    with pytest.raises(SA.QualificationError):
        SA.Surrogate(_ems(), "test", benchmark=bad)                  # benchmark of different emulators
    good = dict(passed=True, emulator_tags=SA._tags(_ems()))
    assert SA.Surrogate(_ems(), "test", benchmark=good)(3.0, 50)["qualified"]
    assert not SA.Surrogate(_ems(), "test", benchmark=dict(good, passed=False), allow_unqualified=True)(3.0, 50)["qualified"]


def test_interpolation_zero_algae_reference_and_flags():
    s = SA.Surrogate(_ems(), "test", allow_unqualified=True)
    r = s(3.25, 50.0, r_um=2500.0, f_n=0.5, dust_ppb=1e5)
    assert r["in_domain"] and not r["flags"]
    assert np.isclose(r["bba"], 0.65 - 0.03 * 3.25 - 0.005 - 1e-6 * 2500 - 0.05)
    # d_alpha from the forcing relative to the algae-free run, per node irradiance: 0.03 x log B exactly
    assert np.isclose(r["dalpha_algae"], 0.03 * 3.25)
    out = s(7.5, 70.0)
    assert not out["in_domain"] and len(out["flags"]) == 2          # log_b and sza; default dust is in range
    assert np.isclose(out["rf_algae"], (800 - 60) * 0.03 * 6.0)


def test_domain_design_covers_ends_and_off_node_sza():
    s = SA.Surrogate(_ems(), "test", allow_unqualified=True)
    st = SA.domain_design(s.ranges, [40.0, 60.0], n_per_sza=8, radii=np.array([1000.0, 4000.0]))
    lbs = {x["log_b"] for x in st}
    assert 1.0 in lbs and 6.0 in lbs and {x["sza"] for x in st} == {40.0, 50.0, 60.0}
    assert {x["r_um"] for x in st} <= {1000.0, 4000.0}
