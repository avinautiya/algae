"""Verification of the independent resolved reference solver itself (it is only usable as a reference if these pass)."""
import os, sys
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "ice_network", "tests"))
from ice_network_reference.solver import CP_W, G, LF, RHO_I, RHO_W, RefNetwork, htc
import reference as R          # depth-ODE steady reference (independent of both solvers)

zero = lambda t: np.array([0.0])


def channel(n, z0=100.4985, drop=1.0, w=0.5, L=100.0, pool=95.0, ke=0.5, nm=0.04):
    zc = z0 - drop * (np.arange(n) + 0.5) / n
    return dict(lake=0, z=zc, w=np.full(n, w), zs=zc + 1.0, length=L, n_man=nm, ke=ke, pool=pool)


def test_lake_at_rest_is_preserved_over_variable_bed_and_width():
    n = 120; x = (np.arange(n) + 0.5) / n
    z = 99.0 + 0.3 * np.sin(6 * x) - 0.2 * x; w = 0.5 + 0.2 * np.sin(9 * x)
    c = dict(lake=0, z=z, w=w, zs=z + 3.0, length=60.0, n_man=0.04, ke=0.5, pool=101.0)
    net = RefNetwork([dict(z_bottom=95.0, area=1e6)], [c], melting=False)
    v0 = 1e6 * (101.0 - 95.0)
    net.run([v0], [0.0], 300.0, zero, zero, record_every=100.0)
    # start from the rest state: depth = 101 - z everywhere
    net.ch[0]["h"][:] = 101.0 - z
    net.ch[0]["Q"][:] = 0.0
    net.run([v0], [0.0], 300.0, zero, zero, record_every=100.0)
    assert np.abs(net.ch[0]["Q"]).max() < 1e-12
    assert np.abs(net.ch[0]["z"] + net.ch[0]["h"] - 101.0).max() < 1e-10


def test_steady_discharge_converges_to_the_independent_depth_ode_solution():
    Hl = 100.72
    with np.errstate(all="ignore"):
        q_ode = R.gvf_ode_Q(Hl, 95.0, 0.5, 100.4985, 99.4985, 100.0, 0.04, 0.5)
    errs = []
    for n in (50, 100, 200):
        net = RefNetwork([dict(z_bottom=90.0, area=1e7)], [channel(n)], melting=False)
        net.run([1e7 * (Hl - 90.0)], [0.0], 2500.0, zero, zero, record_every=2500.0)
        errs.append(abs(net.ch[0]["Q"][n // 2] - q_ode) / q_ode)
    assert errs[0] > errs[1] > errs[2]
    order = np.log2(errs[1] / errs[2])
    assert 0.6 < order < 1.5                      # first-order scheme
    rich = errs[2] - (errs[1] - errs[2])          # Richardson step for a first-order method
    assert abs(rich) < 4e-3


def test_water_mass_is_conserved_with_melt_water_added_and_a_sill_overtopping():
    n = 100
    net = RefNetwork([dict(z_bottom=99.0, area=500.0)], [channel(n, z0=100.0, drop=0.5)], melting=True)
    r = net.run([500.0 * 0.95], [0.4], 3 * 3600.0, lambda t: np.array([0.02]), lambda t: np.array([0.4]), record_every=1800.0)
    last = r["rec"][-1]; first = r["rec"][0]
    stored = last["V"].sum() + last["channel_volume"].sum()
    resid = stored + last["delivered"].sum() - (first["V"].sum() + first["channel_volume"].sum()) - last["input"] - last["melt_w"]
    assert last["delivered"].sum() > 0 and last["melt_w"] > 0
    assert abs(resid) / last["input"] < 1e-9


def test_thermal_energy_closes_in_the_reference_with_explicit_dissipation():
    n = 100
    net = RefNetwork([dict(z_bottom=99.0, area=500.0)], [channel(n, z0=100.0, drop=0.5)], melting=True)
    r = net.run([500.0 * 1.05], [0.4], 2 * 3600.0, lambda t: np.array([0.02]), lambda t: np.array([0.4]), record_every=3600.0)
    f, l = r["rec"][0], r["rec"][-1]
    stored = (l["lake_sens"] + l["chan_sens"]) - (f["lake_sens"] + f["chan_sens"])
    lhs = stored + l["sens_out"] + l["heat"]
    rhs = l["sens_in_ext"] + l["diss"]
    assert l["heat"] > 0 and l["sens_out"] > 0
    assert abs(lhs - rhs) / abs(rhs) < 1e-8


def test_thermal_front_and_cooling_follow_the_independent_plug_flow_solution():
    """Codex 'thermal_front' limit: uniform steady flow, step in inlet temperature, linear wall cooling, T = Tin exp(-beta x/u) behind the front."""
    n = 400; L = 100.0; w = 0.5; nm = 0.04; s0 = 0.005
    hn = 0.2; a = w * hn; rr = a / (w + 2 * hn); qn = a * rr ** (2 / 3) * np.sqrt(s0) / nm; vn = qn / a
    zc = 100.0 - s0 * L * (np.arange(n) + 0.5) / n
    c = dict(lake=0, z=zc, w=np.full(n, w), zs=zc + 1.0, length=L, n_man=nm, ke=0.0, pool=zc[-1] - 0.0 + hn)
    net = RefNetwork([dict(z_bottom=90.0, area=1e8)], [c], melting=True, heat_diss=False)
    net.ch[0]["h"][:] = hn; net.ch[0]["Q"][:] = qn
    h_lake = zc[0] + hn + vn * vn / (2 * G)
    # keep the lake stage fixed by a huge basin; start the lake at the right stage
    v0 = 1e8 * (h_lake - 90.0)
    tin = 1.0
    net.run([v0], [tin], 0.5 * L / vn, zero, zero, record_every=1e9)              # front has travelled half the channel
    beta = htc(qn, w, hn, 1.0) * (w + 2 * hn) / (RHO_W * CP_W * a)
    x = (np.arange(n) + 0.5) * (L / n); T = net.ch[0]["T"]
    front = vn * 0.5 * L / vn                                                      # = L/2
    behind = x < front - 10.0
    exact = tin * np.exp(-beta * x[behind] / vn)
    assert np.abs(T[behind] - exact).max() < 0.03 * tin                          # first-order upwind: small smearing behind the front
    assert T[x > front + 10.0].max() < 5e-3 * tin                                 # only numerical diffusion of the first-order upwind scheme ahead of the front
    # position of the 50 % level of the front is within a few cells of x = u t
    half = 0.5 * tin * np.exp(-beta * front / vn)
    xfront = x[np.argmax(T < half)]
    assert abs(xfront - front) < 2.0


def _two_channel_net():
    n1, n2 = 120, 80
    z1 = 100.0 - 0.01 * 60.0 * (np.arange(n1) + 0.5) / n1; z2 = 100.04 - 0.02 * 40.0 * (np.arange(n2) + 0.5) / n2
    c1 = dict(lake=0, z=z1, w=np.full(n1, 0.4), zs=z1 + 1.0, length=60.0, n_man=0.04, ke=0.5, pool=90.0)
    c2 = dict(lake=0, z=z2, w=np.full(n2, 0.4), zs=z2 + 1.0, length=40.0, n_man=0.04, ke=0.5, pool=90.0)
    return RefNetwork([dict(z_bottom=99.0, area=300.0)], [c1, c2], melting=True)


def test_budgets_close_through_wetting_competing_outlets_drying_and_melting():
    """Two dry channels fed by one lake: the lake rises, wets R1 then R2 (junction with two outlets), melts, then the inflow stops and the
    lake drains below both sills so the channels dry. Water and thermal-energy ledgers must close at EVERY record, not only at the end."""
    net = _two_channel_net()
    inflow = lambda t: np.array([0.03 if t < 1500.0 else 0.0])
    r = net.run([300.0 * 0.97], [0.8], 7200.0, inflow, lambda t: np.array([0.8]), record_every=60.0)
    rec = r["rec"]; f = rec[0]
    wres, eres = [], []
    for x in rec[1:]:
        stored = x["V"].sum() + x["channel_volume"].sum() - f["V"].sum() - f["channel_volume"].sum()
        wres.append(abs(stored + x["delivered"].sum() - x["input"] - x["melt_w"]) / max(x["input"], 1e-9))
        lhs = (x["lake_sens"] + x["chan_sens"]) - (f["lake_sens"] + f["chan_sens"]) + x["sens_out"] + x["heat"]
        rhs = x["sens_in_ext"] + x["diss"]
        if abs(rhs) > 1e6:
            eres.append(abs(lhs - rhs) / abs(rhs))
    q1 = np.array([x["Q_in"][0] for x in rec]); q2 = np.array([x["Q_in"][1] for x in rec])
    assert net.min_depth_seen >= 0.0                                   # positivity through wetting and drying
    assert q1.max() > 1e-3 and q2.max() > 1e-4                         # both routes were wetted (R2 after R1: competing outlets)
    assert np.argmax(q2 > 1e-4) > np.argmax(q1 > 1e-4)
    cv2 = np.array([x["channel_volume"][1] for x in rec])
    assert q2[-1] < 1e-3 * q2.max() and cv2[-1] < 0.01 * cv2.max()      # R2 wetted and then dried while the lake falls below its sill
    assert rec[-1]["stage"][0] < rec[-1]["sill"][1] and q1[-1] < 0.2 * q1.max()   # R1 draining down
    assert rec[-1]["heat"] > 0 and rec[-1]["melt_w"] > 0               # melting happened
    assert max(wres) < 1e-9, max(wres)
    assert max(eres) < 1e-8, max(eres)
    print("max water residual %.2e, max energy residual %.2e, min depth %.2e" % (max(wres), max(eres), net.min_depth_seen))
