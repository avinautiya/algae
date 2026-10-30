"""Fixed-geometry hydraulics: analytical and independent-numerical references, refinement, backwater, reversal."""
import numpy as np
import pytest

from ice_network.hydraulics.channel import edge_flow
from ice_network.constants import G
import reference as R

NM, KE = 0.04, 0.5


def call(h_i, h_j, z, w, dx, msub=16, nm=NM, ke=KE):
    n = z.size
    h = np.zeros(n); hf = np.zeros(n - 1); ch = np.zeros(n, np.int64)
    out = edge_flow(h_i, h_j, z, w, nm, dx, ke, msub, h, hf, ch)
    return out, h, hf, ch


def channel(z0, drop, w0, L=100.0, n=9):
    return z0 - drop * np.linspace(0, 1, n), np.full(n, w0), L / (n - 1)


def test_uniform_flow_recovers_manning_normal_discharge():
    z, w, dx = channel(100.0, 0.5 * 4.0 / 4.0 * 1.0 * 0 + 2.0, 0.6, L=400.0, n=21)   # S0 = 0.005
    s0 = 2.0 / 400.0
    hn = 0.2
    a = 0.6 * hn; r = a / (0.6 + 2 * hn); qn = a * r ** (2 / 3) * np.sqrt(s0) / NM
    v = qn / a
    h_res = z[0] + hn + (1 + KE) * v * v / (2 * G)
    (q, d, _, _), h, _, ch = call(h_res, z[-1] + hn, z, w, dx)
    assert q == pytest.approx(qn, rel=1e-9)
    assert np.ptp(h) < 1e-9 and ch.sum() == 0 and d == 1


@pytest.mark.parametrize("case", [
    dict(name="mild free outfall", z0=100.4985, drop=1.0, hres=100.72, tail=95.0, tol=1e-6),
    dict(name="mild submerged", z0=100.4985, drop=1.0, hres=100.80, tail=99.60, tol=1e-6),
    dict(name="flat ponded near-critical", z0=100.0, drop=0.02, hres=100.30, tail=99.99, tol=5e-3),
])
def test_against_independent_depth_ode_reference(case):
    z, w, dx = channel(case["z0"], case["drop"], 0.5)
    with np.errstate(all="ignore"):
        q_ref = R.gvf_ode_Q(case["hres"], case["tail"], 0.5, z[0], z[-1], 100.0, NM, KE)
    (q, *_), *_ = call(case["hres"], case["tail"], z, w, dx)
    assert abs(q - q_ref) / q_ref < case["tol"]


def test_substep_refinement_converges_at_second_order():
    z, w, dx = channel(100.4985, 1.0, 0.5)
    with np.errstate(all="ignore"):
        q_ref = R.gvf_ode_Q(100.80, 99.60, 0.5, z[0], z[-1], 100.0, NM, KE)
    err = [abs(call(100.80, 99.60, z, w, dx, msub=m)[0][0] - q_ref) / q_ref for m in (2, 4, 8, 16)]
    assert all(e2 < e1 for e1, e2 in zip(err, err[1:]))
    orders = [np.log2(e1 / e2) for e1, e2 in zip(err, err[1:])]
    assert all(1.7 < o < 2.3 for o in orders), orders


def test_inlet_control_matches_critical_flow_weir_law():
    z = np.array([101.0, 100.0]); w = np.array([0.5, 0.5])
    for hres in (101.1, 101.3, 101.6):
        (q, *_), *_ = call(hres, 95.0, z, w, 10.0, nm=1e-4)
        hc = (hres - 101.0) / (1.5 + 0.5 * KE)
        assert q == pytest.approx(0.5 * np.sqrt(G) * hc ** 1.5, rel=1e-9)


def test_backwater_reduces_discharge_and_matches_reference_on_a_short_channel():
    # outlet-control regime: modest head over the sill, small bed drop, so friction and tailwater govern
    z, w, dx = channel(100.0, 0.05, 0.5, L=20.0, n=5)
    tails = 99.95 + np.array([0.0, 0.05, 0.10, 0.15, 0.20, 0.24])
    qs = []
    for t in tails:
        with np.errstate(all="ignore"):
            q_ref = R.gvf_ode_Q(100.20, t, 0.5, z[0], z[-1], 20.0, NM, KE)
        (q, *_), *_ = call(100.20, t, z, w, dx, msub=32)
        qs.append(q)
        assert abs(q - q_ref) / q_ref < 2e-3
    # free outfall (tailwater depth below critical) leaves Q unchanged; once the exit is submerged Q must fall
    assert all(a >= b for a, b in zip(qs, qs[1:])), qs
    assert qs[3] > qs[4] > qs[5]
    assert qs[0] / qs[-1] > 1.01


def test_inlet_controlled_discharge_is_insensitive_to_tailwater_until_drowned():
    z, w, dx = channel(100.4985, 1.0, 0.5)
    q = [call(100.72, t, z, w, dx)[0][0] for t in (95.0, 97.0, 98.0)]
    assert max(q) / min(q) - 1 < 1e-9


def test_direction_follows_the_water_surface_and_is_antisymmetric():
    z, w, dx = channel(100.0, 0.2, 0.5, L=50.0, n=6)
    (qf, df, *_), *_ = call(100.5, 100.1, z, w, dx)
    (qr, dr, *_), *_ = call(100.1, 100.5, z[::-1].copy() * 0 + z, w, dx)
    assert df == 1 and qf > 0 and dr == -1 and qr < 0
    (q0, *_), *_ = call(100.3, 100.3, z, w, dx)
    assert q0 == 0.0


def test_dry_channel_and_hump_block_flow():
    z, w, dx = channel(100.0, 0.2, 0.5, L=50.0, n=6)
    assert call(99.9, 95.0, z, w, dx)[0][0] == 0.0               # reservoir below the sill
    zh = z.copy(); zh[3] = 100.6                                  # an ice hump higher than the reservoir stage
    assert call(100.5, 95.0, zh, w, dx)[0][0] == 0.0
    assert call(100.7, 95.0, zh, w, dx)[0][0] > 0.0
