"""Thermal transport and melting: analytical references, energy/latent-heat balance, kilogram accounting."""
import numpy as np
import pytest
from scipy.optimize import brentq

from ice_network.constants import CP_W, K_W, LF, NU_W, PR_W, RHO_I, RHO_W
from ice_network.evolution.network import Network
from ice_network.thermal.channel_heat import reach_heat
from helpers import chan, const, lake, q_of_stage

RC = RHO_W * CP_W


def htc_independent(q, w, h):
    a, p = w * h, w + 2 * h
    dh = 4 * a / p
    re = (q / a) * dh / NU_W
    return 0.023 * re ** 0.8 * PR_W ** 0.3 * K_W / dh


@pytest.mark.parametrize("n", [2, 5, 17])
def test_reach_temperature_decays_exponentially_independent_of_station_count(n):
    q, w0, h0, L, t_in = 0.05, 0.5, 0.2, 100.0, 0.3
    h = np.full(n, h0); w = np.full(n, w0); dx = L / (n - 1)
    t_arr = np.zeros(n); qi = np.zeros(n - 1)
    t_out = reach_heat(q, t_in, 0.0, h, w, dx, np.zeros(n - 1), 0.0, 1.0, t_arr, qi)
    htp = htc_independent(q, w0, h0) * (w0 + 2 * h0)
    t_exact = t_in * np.exp(-htp * L / (RC * q))
    assert t_out == pytest.approx(t_exact, rel=1e-12)
    assert qi.sum() == pytest.approx(RC * q * (t_in - t_exact), rel=1e-11)         # heat to ice = sensible loss


def test_reach_with_dissipation_source_matches_closed_form_and_heat_integral():
    q, w0, h0, L, t_in, n = 0.05, 0.5, 0.2, 100.0, 0.05, 9
    dx = L / (n - 1); phi_line = 40.0                                                # W per metre of channel
    h = np.full(n, h0); w = np.full(n, w0)
    t_arr = np.zeros(n); qi = np.zeros(n - 1)
    t_out = reach_heat(q, t_in, 0.0, h, w, dx, np.full(n - 1, phi_line * dx), 0.0, 1.0, t_arr, qi)
    htp = htc_independent(q, w0, h0) * (w0 + 2 * h0)
    lam = RC * q / htp; t_eq = phi_line / htp
    assert t_out == pytest.approx(t_eq + (t_in - t_eq) * np.exp(-L / lam), rel=1e-12)
    heat_exact = htp * (t_eq * L + (t_in - t_eq) * lam * (1 - np.exp(-L / lam)))      # h_T P integral of T dx
    assert qi.sum() == pytest.approx(heat_exact, rel=1e-11)


def test_melt_rate_limit_and_geometry_update_conserve_ice_volume_on_a_uniform_segment():
    q, w0, h0, dx, t_in = 5.0, 0.5, 0.2, 10.0, 0.3        # large Q: negligible cooling over the segment
    h = np.full(2, h0); w = np.full(2, w0); t_arr = np.zeros(2); qi = np.zeros(1)
    reach_heat(q, t_in, 0.0, h, w, dx, np.zeros(1), 0.0, 1.0, t_arr, qi)
    p = w0 + 2 * h0
    v_melt = qi[0] / (p * dx * RHO_I * LF)
    htc = htc_independent(q, w0, h0)
    v_stefan = htc * t_in / (RHO_I * LF)                                                # h_T (T_in - T_m) / (rho_i L_f)
    x = htc * p * dx / (RC * q)                                                          # cooling over the segment
    assert x < 0.02 and v_melt < v_stefan                                                # the inlet-temperature limit is an upper bound
    assert v_melt == pytest.approx(v_stefan * (1 - np.exp(-x)) / x, rel=1e-10)         # exact segment-mean temperature
    e = chan(0, 1, 100.0, 0.1, w0, L=dx, n=2)
    net = Network([lake()], [dict(h_fixed=95.0)], [e])
    dt = 3600.0
    wall = 2 * h0 / p
    ve, vg = net.update_geometry(np.array([qi[0] * dt]), np.array([qi[0] * dt * wall]))
    assert ve == pytest.approx(qi[0] * dt / (RHO_I * LF), rel=1e-14)
    assert vg == pytest.approx(ve, rel=1e-10)                                          # geometry change = energy-based volume
    assert net.edges[0]["z"][0] < 100.0 and net.edges[0]["w"][0] > w0


def _steady_two_tank(I):
    e1 = chan(0, 1, 100.0, 0.3, 0.5, L=50.0, n=6, msub=8); e2 = chan(1, 2, 99.0, 0.3, 0.5, L=50.0, n=6, msub=8)
    h2 = brentq(lambda H: q_of_stage(e2, H, 95.0) - I, 99.0 + 1e-9, 100.5, xtol=1e-14)
    h1 = brentq(lambda H: q_of_stage(e1, H, h2) - I, 100.0 + 1e-9, 101.5, xtol=1e-14)
    return e1, e2, h1, h2


def test_thermal_step_advects_and_mixes_through_two_tanks_like_the_analytic_chain():
    I, T0 = 0.05, 0.4
    e1, e2, h1, h2 = _steady_two_tank(I)
    A1, A2 = 400.0, 250.0
    lakes = [lake(98.0, A1), lake(97.5, A2)]
    net = Network(lakes, [dict(h_fixed=95.0)], [e1, e2], melting=False, heat_dissipation=False)
    V1, V2 = A1 * (h1 - 98.0), A2 * (h2 - 97.5)
    r = net.run([V1, V2], [0.0, 0.0], 6 * 3600.0, 900.0, const([I, 0.0]), const([T0, 0.0]), evolve=False, rtol=1e-11, atol=1e-12)
    t1, t2 = V1 / I, V2 / I
    for t, s in zip(r["t"], r["hist"]):
        a1 = T0 * (1 - np.exp(-t / t1))
        a2 = T0 * (1 - (t1 * np.exp(-t / t1) - t2 * np.exp(-t / t2)) / (t1 - t2))
        assert s["T"][0] == pytest.approx(a1, abs=2e-6 * T0)
        assert s["T"][1] == pytest.approx(a2, abs=2e-6 * T0)
    b = r["budgets"]
    assert abs(b["energy_residual_rel"]) < 1e-9 and abs(b["mass_residual_rel"]) < 1e-9
    assert b["heat_to_ice_J"] == 0.0


def _melting_case(n=9, dt=900.0, rtol=1e-9, t_end=6 * 3600.0, melting=True, evolve=True, I=0.05, T0=0.3):
    e = chan(0, 1, 100.5, 1.0, 0.5, L=100.0, n=n)
    net = Network([lake(99.0, 2000.0)], [dict(h_fixed=95.0)], [e], melting=melting)
    r = net.run([2000.0 * 1.55], [T0], t_end, dt, const(I), const(T0), rtol=rtol, atol=1e-10, evolve=evolve)
    return net, r


def test_melting_case_latent_heat_matches_available_energy_and_mass_in_kg():
    net, r = _melting_case()
    b = r["budgets"]
    assert b["delivered_m3"][0] > 100.0 and b["heat_to_ice_J"] > 0                           # nonzero throughput, nonzero melt
    assert abs(b["water_residual_rel"]) < 1e-8 and abs(b["mass_residual_rel"]) < 1e-8
    assert abs(b["energy_residual_rel"]) < 1e-8
    # three independent routes to the melted ice mass (kg):
    m_heat, m_water, m_geom = b["ice_mass_from_heat_kg"], b["ice_mass_from_meltwater_kg"], b["ice_mass_from_geometry_kg"]
    assert m_water == pytest.approx(m_heat, rel=1e-9)            # meltwater added to the flow = latent-heat mass
    assert m_geom == pytest.approx(m_heat, rel=2e-3)             # geometry change vs energy: O(dx^2, dt) agreement (see convergence study)
    latent = RHO_I * LF * b["ice_volume_geometry"]
    assert latent == pytest.approx(b["heat_to_ice_J"], rel=2e-3)


def test_no_melting_control_keeps_geometry_and_still_conserves():
    net, r = _melting_case(melting=False, evolve=False)
    b = r["budgets"]
    assert b["heat_to_ice_J"] == 0.0 and b["melt_water_m3"] == 0.0
    e = net.edges[0]
    assert e["z"][0] == pytest.approx(100.5, abs=0) and np.allclose(e["w"], 0.5)
    assert abs(b["water_residual_rel"]) < 1e-9 and abs(b["energy_residual_rel"]) < 1e-9 and b["dissipation_J"] > 0
