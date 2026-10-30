"""The fast discharge solver (safeguarded Newton stations + Brent on the entrance equation) against the unchanged bisection
reference; residuals, iteration counts, near-dry and reversal behaviour; dissipative-heating switch."""
import numpy as np
import pytest

from ice_network.constants import CP_W, G, RHO_W
from ice_network.evolution.network import Network
from ice_network.hydraulics.channel import edge_flow_diag
from ice_network.hydraulics.channel_ref import edge_flow_ref
from helpers import chan, const, lake

NM, KE = 0.04, 0.5
Q_AMBIGUOUS = 1e-8          # m3/s


def both(h_i, h_j, z, w, dx, msub=16):
    n = z.size
    h = np.zeros(n); hf = np.zeros(n - 1); ch = np.zeros(n, np.int64); dg = np.zeros(4)
    qf = edge_flow_diag(h_i, h_j, z, w, NM, dx, KE, msub, h, hf, ch, dg)[0]
    qr = edge_flow_ref(h_i, h_j, z, w, NM, dx, KE, msub, h, hf, ch)[0]
    return qf, qr, dg.copy()


def test_fast_solver_reproduces_the_reference_over_levels_tailwaters_and_geometries():
    rows = []
    for (z0, drop, w0, L, n) in [(100.4985, 1.0, 0.5, 100.0, 9), (100.0, 0.02, 0.5, 100.0, 9), (100.0, 0.05, 0.3, 20.0, 5)]:
        z = z0 - drop * np.linspace(0, 1, n); w = np.full(n, w0); dx = L / (n - 1)
        for head in np.concatenate([[1e-8, 1e-6, 1e-4], np.geomspace(1e-3, 2.0, 12)]):
            for tail_depth in (-5.0, 0.0, 0.3 * head, 0.8 * head):
                h_i = z0 + head
                h_j = z[-1] + tail_depth
                if h_j >= h_i:
                    continue
                qf, qr, dg = both(h_i, h_j, z, w, dx)
                rows.append((head, qf, qr, dg))
                if qr > Q_AMBIGUOUS:
                    assert abs(qf - qr) / qr < 1e-8, (head, tail_depth, qf, qr)
                else:
                    # below ~1e-8 m3/s the entrance equation is not monotone in Q (micrometric critical depth, Manning slope >> 1), the
                    # root is not unique and the two solvers may select different roots; the absolute difference is irrelevant (<1e-9 m3/s)
                    assert abs(qf - qr) < 1e-9
    rel = [abs(q - r) / r for _, q, r, _ in rows if r > Q_AMBIGUOUS]
    its = [d[0] for _, q, r, d in rows if d[3] == 2.0]
    res = [d[1] for _, q, r, d in rows if d[3] == 2.0]
    assert max(its) <= 30 and np.mean(its) < 12          # observed: mean ~8, max 22 (Brent with a bisection safeguard)
    assert max(res) < 1e-12                  # residual of the entrance equation, metres of stage
    assert max(rel) < 1e-8


def test_near_dry_discharge_vanishes_continuously_and_monotonically():
    z = 100.0 - 0.05 * np.linspace(0, 1, 5); w = np.full(5, 0.5)
    heads = np.geomspace(1e-9, 1e-2, 40)
    qs = np.array([both(100.0 + h, 95.0, z, w, 5.0)[0] for h in heads])
    assert np.all(qs >= 0) and np.all(np.diff(qs) >= 0)
    assert qs[0] < 1e-12 and qs[-1] > 1e-6
    m = qs > 1e-10                              # at micro-heads the friction-limited channel carries no flow at all (Q = 0 is reported)
    assert m.sum() >= 5 and not m[0]
    slope = np.diff(np.log(qs[m])) / np.diff(np.log(heads[m]))
    assert np.all(slope > 0.9) and np.all(slope < 2.2)      # smooth power-law-like onset; no jumps


def test_flow_reversal_is_antisymmetric_and_continuous_through_zero_head_difference():
    n = 5; z = np.full(n, 99.8); w = np.full(n, 0.5); dx = 7.5
    base = 100.1
    diffs = np.linspace(-0.2, 0.2, 21)
    q = []
    for d in diffs:
        h = np.zeros(n); hf = np.zeros(n - 1); ch = np.zeros(n, np.int64)
        from ice_network.hydraulics.channel import edge_flow
        q.append(edge_flow(base - d / 2, base + d / 2, z, w, NM, dx, KE, 16, h, hf, ch)[0])
    q = np.array(q)
    assert q[10] == 0.0
    assert np.allclose(q, -q[::-1], rtol=1e-9, atol=1e-13)          # antisymmetry (flat symmetric channel)
    assert np.all(np.diff(q) < 0)                                     # monotone through zero: positive when a is higher
    sl = np.log(abs(q[0]) / abs(q[9])) / np.log(0.2 / 0.02)
    assert 0.4 < sl < 0.6                                              # Q ~ |dH|^(1/2): friction-controlled, continuous through zero


def test_symmetric_branches_receive_identical_discharge_at_every_level():
    n = 7; z = 100.0 - 0.3 * np.linspace(0, 1, n); w = np.full(n, 0.4)
    for head in np.geomspace(1e-4, 1.0, 10):
        a = both(100.0 + head, 95.0, z, w, 10.0)[0]
        b = both(100.0 + head, 95.0, z.copy(), w.copy(), 10.0)[0]
        assert a == b


# ----------------------------------------------------------------------------------------- dissipative heating switch
def _flow_case(heat_dissipation, melting=False, evolve=False, T0=0.0, t_end=4 * 3600.0, dt=900.0):
    e = chan(0, 1, 100.0, 0.6, 0.4, L=80.0, n=7, msub=8)
    net = Network([lake(99.0, 600.0)], [dict(h_fixed=95.0)], [e], melting=melting, heat_dissipation=heat_dissipation)
    r = net.run([600.0 * 1.4], [T0], t_end, dt, const(0.05), const(T0), evolve=evolve, rtol=1e-10, atol=1e-11)
    return net, r


def test_dissipation_switch_does_not_change_hydraulics_and_adds_exactly_the_dissipated_energy():
    n_on, r_on = _flow_case(True)
    n_off, r_off = _flow_case(False)
    for a, b in zip(r_on["hist"], r_off["hist"]):
        assert np.allclose(a["V"], b["V"], rtol=1e-8, atol=1e-9) and np.allclose(a["Q"], b["Q"], rtol=1e-7, atol=1e-12)   # same friction law, same discharge/storage up to integrator step selection (the extra temperature states change the step sequence)
    bo, bf = r_on["budgets"], r_off["budgets"]
    assert bf["dissipation_J"] == 0.0 and bo["dissipation_J"] > 0.0
    rc = RHO_W * CP_W
    sens_on = rc * 0 + bo["sensible_out_J"] + bo["sensible_storage_change_J"]
    sens_off = bf["sensible_out_J"] + bf["sensible_storage_change_J"]
    assert (sens_on - sens_off) == pytest.approx(bo["dissipation_J"], rel=1e-8)          # the whole difference is the dissipation
    for b in (bo, bf):
        assert abs(b["energy_residual_rel"]) < 1e-9 and abs(b["water_residual_rel"]) < 1e-9


def test_steady_outflow_temperature_rise_equals_g_times_total_head_drop():
    net, r = _flow_case(True, t_end=12 * 3600.0)
    s = r["hist"][-1]
    stage_drop = s["stage"][0] - s["stage"][1]
    q = abs(s["Q"][0])
    # late-time sensible heat leaving per second (finite difference of the cumulative outlet heat)
    n2, r2 = _flow_case(True, t_end=12 * 3600.0 + 900.0)
    de = (r2["budgets"]["sensible_out_J"] - r["budgets"]["sensible_out_J"]) / 900.0
    assert de == pytest.approx(RHO_W * G * q * stage_drop, rel=2e-3)                    # T_out - T_in = g dH / c_w at steady state


def test_dissipation_only_melting_never_exceeds_the_channel_share_of_dissipation():
    net, r = _flow_case(True, melting=True, evolve=True, T0=0.0, t_end=6 * 3600.0)
    b = r["budgets"]
    assert 0.0 < b["heat_to_ice_J"] <= b["dissipation_channel_J"] * (1 + 1e-9)
    assert b["dissipation_exit_J"] > 0.0                                                # exit losses heat the outlet water, not the ice
    assert abs(b["energy_residual_rel"]) < 1e-8
