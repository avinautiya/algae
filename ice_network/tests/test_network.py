"""Storage, branching, symmetry, activation, reversal; mass (m3 and kg) conservation on nonzero-throughput cases."""
import numpy as np
import pytest
from scipy.optimize import brentq

from ice_network.constants import G, RHO_W
from ice_network.evolution.network import Network
from helpers import chan, const, lake, q_of_stage


def run(net, v0, t0, t_end, dt, inflow, tin, **kw):
    return net.run(v0, t0, t_end, dt, inflow, tin, evolve=False, **kw)


def test_reservoir_drains_through_inlet_controlled_channel_like_the_analytic_solution():
    A, ke, w = 1000.0, 0.5, 0.5
    e = chan(0, 1, 100.0, 1.0, w, L=10.0, n=2, nm=1e-4, ke=ke)
    net = Network([lake(99.0, A)], [dict(h_fixed=95.0)], [e], melting=False)
    h0 = 0.5
    r = run(net, [A * (1.0 + h0)], [0.0], 3600.0, 600.0, const(0.0), const(0.0), rtol=1e-11, atol=1e-12)
    c = w * np.sqrt(G) / (A * (1.5 + 0.5 * ke) ** 1.5)
    for t, snap in zip(r["t"], r["hist"]):
        h_exact = (h0 ** -0.5 + 0.5 * c * t) ** -2
        assert snap["stage"][0] - 100.0 == pytest.approx(h_exact, rel=2e-6)
    b = r["budgets"]
    assert abs(b["water_residual_rel"]) < 1e-9 and abs(b["mass_residual_rel"]) < 1e-9
    assert b["delivered_m3"][0] > 0.2 * A * h0       # nonzero throughput


def test_branching_junction_steady_split_and_budgets():
    I = 0.06
    e1 = chan(0, 1, 100.0, 0.5, 0.3, L=60.0, n=7)
    e2 = chan(0, 2, 100.05, 0.8, 0.5, L=80.0, n=9)
    net = Network([lake(99.0, 500.0)], [dict(h_fixed=95.0), dict(h_fixed=94.0)], [e1, e2], melting=False)
    f = lambda H: q_of_stage(e1, H, 95.0) + q_of_stage(e2, H, 94.0) - I
    Hs = brentq(f, 100.06, 101.5, xtol=1e-14)
    r = run(net, [500.0 * (Hs + 0.15 - 99.0)], [0.0], 24 * 3600.0, 1800.0, const(I), const(0.0))
    last = r["hist"][-1]
    assert last["stage"][0] == pytest.approx(Hs, abs=1e-6)
    assert last["Q"][0] == pytest.approx(q_of_stage(e1, Hs, 95.0), rel=1e-4)
    assert last["Q"][1] == pytest.approx(q_of_stage(e2, Hs, 94.0), rel=1e-4)
    assert last["Q"].sum() == pytest.approx(I, rel=1e-4)
    b = r["budgets"]
    assert abs(b["water_residual_rel"]) < 1e-9 and abs(b["mass_residual_rel"]) < 1e-9
    assert b["delivered_m3"].sum() > 0.8 * b["input_m3"]          # nonzero throughput: most input leaves
    assert b["mass_out_kg"] == pytest.approx(RHO_W * b["delivered_m3"].sum(), rel=1e-12)


def test_symmetric_outlets_receive_equal_discharge():
    e1 = chan(0, 1, 100.0, 0.5, 0.4); e2 = chan(0, 2, 100.0, 0.5, 0.4)
    net = Network([lake(99.0, 500.0)], [dict(h_fixed=95.0), dict(h_fixed=95.0)], [e1, e2], melting=False)
    r = run(net, [500.0 * 1.3], [0.0], 6 * 3600.0, 1200.0, const(0.08), const(0.0))
    c = r["budgets"]["delivered_m3"]
    assert c[0] > 100.0 and c[0] == pytest.approx(c[1], rel=1e-9)


def test_dry_to_wet_activation_keeps_storage_nonnegative():
    e = chan(0, 1, 100.0, 0.5, 0.4)
    net = Network([lake(99.0, 500.0)], [dict(h_fixed=95.0)], [e], melting=False)
    r = run(net, [1e-3], [0.0], 6 * 3600.0, 300.0, const(0.05), const(0.0))
    sill = 100.0
    assert min(s["V"][0] for s in r["hist"]) >= 0.0
    first_flow = next(i for i, s in enumerate(r["hist"]) if s["Q"][0] > 0)
    assert first_flow > 5
    assert all(s["stage"][0] <= sill + 1e-9 for s in r["hist"][:first_flow])
    assert all(s["stage"][0] > sill for s in r["hist"][first_flow:])
    assert abs(r["budgets"]["water_residual_rel"]) < 1e-9


def test_direction_reversal_between_two_storage_nodes_and_equalisation():
    e = chan(0, 1, 99.8, 0.0, 0.5, L=30.0, n=5)
    e["z"] = np.full(5, 99.8); e["zs"] = e["z"] + 1.0
    net = Network([lake(99.0, 500.0), lake(99.0, 500.0)], [], [e], melting=False)
    I = lambda t: np.array([0.3 if t < 1800.0 else 0.0, 0.0])
    v0 = [500.0 * 1.2, 500.0 * 1.5]
    r = run(net, v0, [0.0, 0.0], 8 * 3600.0, 300.0, I, const([0.0, 0.0]), rtol=1e-10, atol=1e-11)
    sgn = np.sign([s["Q"][0] for s in r["hist"]])
    first_nonzero = sgn[sgn != 0]
    assert first_nonzero[0] == -1                       # b -> a initially (b higher)
    assert (first_nonzero == 1).any()                   # reverses once a is raised above b
    assert abs(r["hist"][-1]["stage"][0] - r["hist"][-1]["stage"][1]) < 1e-3
    b = r["budgets"]
    assert abs(b["water_residual_rel"]) < 1e-9 and b["delivered_m3"].size == 0
    tot = sum(r["hist"][-1]["V"])
    assert tot == pytest.approx(sum(v0) + 0.3 * 1800.0, rel=1e-9)
