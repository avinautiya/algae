"""Treatments used in the controlled incision comparison: dissipation-all melting and fixed-aspect geometry (explicit approximations)."""
import numpy as np
import pytest

from ice_network.constants import LF, RHO_I
from ice_network.evolution.network import Network
from helpers import chan, const, lake


def case(melt_mode="thermal", geometry_mode="energy_partition", T0=0.0, evolve=True, melting=True):
    e = chan(0, 1, 100.5, 1.0, 0.4, L=100.0, n=9)
    net = Network([lake(99.0, 2000.0)], [dict(h_fixed=95.0)], [e], melting=melting, melt_mode=melt_mode, geometry_mode=geometry_mode)
    r = net.run([2000.0 * 1.55], [T0], 4 * 3600.0, 900.0, const(0.05), const(T0), rtol=1e-10, atol=1e-11, evolve=evolve)
    return net, r


def test_dissipation_all_melts_exactly_the_channel_dissipation_and_closes_the_energy_budget():
    net, r = case("dissipation_all", T0=0.0, evolve=False)
    b = r["budgets"]
    assert b["heat_to_ice_J"] == pytest.approx(b["dissipation_channel_J"], rel=1e-9)
    assert b["dissipation_exit_J"] > 0
    assert abs(b["energy_residual_rel"]) < 1e-9 and abs(b["water_residual_rel"]) < 1e-9


def test_thermal_treatment_melts_less_than_dissipation_all_when_exposure_is_modest():
    _, r_t = case("thermal", T0=0.0, evolve=False)
    _, r_d = case("dissipation_all", T0=0.0, evolve=False)
    assert r_t["budgets"]["heat_to_ice_J"] < r_d["budgets"]["heat_to_ice_J"]


def test_fixed_aspect_geometry_keeps_the_ratio_and_the_ice_volume():
    net, r = case("dissipation_all", "fixed_aspect", T0=0.0)
    e = net.edges[0]
    r0 = 0.4 / 1.0
    ratio = e["w"] / (e["zs"] - e["z"])
    assert np.allclose(ratio, r0, rtol=1e-6)
    b = r["budgets"]
    assert b["ice_volume_geometry"] == pytest.approx(b["ice_volume_energy"], rel=1e-3)
