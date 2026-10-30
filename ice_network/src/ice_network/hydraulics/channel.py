"""Quasi-steady open-channel hydraulics of one reach connecting two storage nodes.

Formulation (see docs/mathematical_contract.md, section 3)
  * Rectangular section of width w_k and bed elevation z_k at stations k = 0..N-1 (spacing dx), Manning friction.
  * Steady gradually varied flow solved by the standard-step method on the TOTAL HEAD E = z + h + Q^2/(2 g w^2 h^2):
        E_k = E_{k+1} + 0.5 (Sf_k + Sf_{k+1}) dx          (marching upstream, subcritical branch)
    from a downstream depth h_L = max(H_j - z_{N-1}, h_c) (tailwater or critical depth at the brink).
  * Control points: where the head cannot supply the critical specific energy the depth is set to critical depth and
    the energy line is raised to z + 1.5 h_c (weir / hump / brink control). The raised amount is NOT counted as friction.
  * Entrance: H_i = E_0 + k_e v_0^2 / (2 g) (reservoir velocity head neglected).
  * Each segment is subdivided into `msub` substeps; convergence in msub is verified (tests/test_hydraulics.py).
  * Capacity: Q = min(Q_outlet, Q_inlet), Q_outlet from the root of H_req(Q) = H_i, Q_inlet from critical flow at
    the entrance, h_c = (H_i - z_0) / (1.5 + 0.5 k_e).
  * Direction follows the water-surface gradient; the reverse direction uses the reversed station arrays.
  * Momentum is not solved (no hydraulic jumps, no waves). Supercritical reaches are NOT represented: they are
    reported as 'controlled' stations (depth set to critical) and flagged.
"""
import numpy as np
from numba import njit

from ..constants import G


@njit(cache=True)
def _crit_depth(q, w):
    return (q * q / (G * w * w)) ** (1.0 / 3.0)


@njit(cache=True)
def _sf(q, w, h, n_man):
    a = w * h
    r = a / (w + 2.0 * h)
    return (n_man * q / (a * r ** (2.0 / 3.0))) ** 2


@njit(cache=True)
def _station_depth(e_next, sf_next, z, w, q, n_man, dx):
    """Implicit standard-step station. Governing equation (subcritical branch, h > h_c):
        g(h) = z + h + q^2/(2 g w^2 h^2) - e_next - 0.5 (Sf(h) + sf_next) dx = 0,
        g'(h) = 1 - q^2/(g w^2 h^3) - 0.5 dx dSf/dh,   dSf/dh = Sf (-(10/3)/h + (8/3)/P).
    Safeguarded Newton: the root stays bracketed in [h_c, h_hi] (g(h_c) < 0 < g(h_hi)); a Newton step that leaves the
    bracket or has g' <= 0 is replaced by bisection. Converged when |dh| <= 1e-14 h (also |g| < 1e-15 (1 + e_next - z)).
    If g(h_c) >= 0 no subcritical depth exists: depth is critical and the energy line is raised to z + 1.5 h_c.
    Returns (h, E_k, choked)."""
    q2 = q * q / (w * w)
    hc = (q2 / G) ** (1.0 / 3.0)
    sf_c = _sf(q, w, hc, n_man)
    g_c = z + hc + q2 / (2.0 * G * hc * hc) - e_next - 0.5 * (sf_c + sf_next) * dx
    if g_c >= 0.0:
        return hc, z + 1.5 * hc, True
    lo = hc
    hi = e_next + 0.5 * (sf_c + sf_next) * dx - z
    if hi <= lo:
        hi = 2.0 * lo
    h = hi
    for _ in range(80):
        sf = _sf(q, w, h, n_man)
        gv = z + h + q2 / (2.0 * G * h * h) - e_next - 0.5 * (sf + sf_next) * dx
        if gv > 0.0:
            hi = h
        else:
            lo = h
        p = w + 2.0 * h
        dsf = sf * (-(10.0 / 3.0) / h + (8.0 / 3.0) / p)
        dg = 1.0 - q2 / (G * h ** 3) - 0.5 * dx * dsf
        hn = h - gv / dg if dg > 0.0 else -1.0
        if not (lo < hn < hi):
            hn = 0.5 * (lo + hi)
        if abs(hn - h) <= 1e-14 * hn:
            h = hn
            break
        h = hn
    return h, e_next + 0.5 * (_sf(q, w, h, n_man) + sf_next) * dx, False


@njit(cache=True)
def gvf_upstream(q, h_l, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out):
    """Standard-step profile for discharge q > 0 given the downstream depth. Each station-to-station segment is
    divided into `msub` equal substeps (geometry linearly interpolated), because the trapezoid friction term is
    inaccurate on a coarse step near critical flow. Returns the reservoir stage that supplies the entrance:
    E_0 + ke v0^2/2g. Fills h_out (depth at the stations), hf_out (friction head loss per segment, zero for
    substeps where the energy line was raised at a control) and choked_out (1 if any substep of the segment, or the
    station itself, was a control)."""
    n = z.size
    q2 = q * q
    h_out[n - 1] = h_l
    choked_out[n - 1] = 0
    e_next = z[n - 1] + h_l + q2 / (2.0 * G * (w[n - 1] * h_l) ** 2)
    sf_next = _sf(q, w[n - 1], h_l, n_man)
    dxs = dx / msub
    h_k = h_l
    for k in range(n - 2, -1, -1):
        hf = 0.0
        ch_any = 0
        for j in range(1, msub + 1):
            f = j / msub
            zz = z[k + 1] + (z[k] - z[k + 1]) * f
            ww = w[k + 1] + (w[k] - w[k + 1]) * f
            h_k, e_k, ch = _station_depth(e_next, sf_next, zz, ww, q, n_man, dxs)
            sf_k = _sf(q, ww, h_k, n_man)
            if ch:
                ch_any = 1
            else:
                hf += 0.5 * (sf_k + sf_next) * dxs
            e_next = e_k
            sf_next = sf_k
        h_out[k] = h_k
        choked_out[k] = ch_any
        hf_out[k] = hf
    v0 = q / (w[0] * h_out[0])
    return e_next + ke * v0 * v0 / (2.0 * G)


@njit(cache=True)
def _req(q, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out):
    n = z.size
    hl = max(h_j - z[n - 1], _crit_depth(q, w[n - 1]))
    return gvf_upstream(q, hl, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out)


@njit(cache=True)
def edge_flow_forward(h_i, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out, diag):
    """Discharge from node i to node j for h_i >= h_j. Returns (Q, E_exit_total, E_entry_total), Q >= 0.

    Governing equation: R(Q) = H_req(Q) - h_i = 0 on [q_lo, q_in], where H_req(Q) = E_0(Q) + k_e v_0^2/(2g) comes from the
    standard-step profile and R is increasing in Q; if R(q_in) <= 0 the entrance (critical-flow) control governs, Q = q_in.
    Solver: Brent's method (bracketing + inverse quadratic/secant steps), converged when the bracket is below
    2e-15 Q or |R| <= 1e-14 (1 + h_i - z_0) m. diag = (iterations, |R| at the returned Q [m], bracket width/Q, 0/1/2 =
    dry/entrance control/outlet control)."""
    n = z.size
    diag[0] = 0.0; diag[1] = 0.0; diag[2] = 0.0; diag[3] = 0.0
    if h_i <= z[0] or h_i <= h_j:
        return 0.0, 0.0, 0.0
    zmax = z[0]
    for k in range(n):
        if z[k] > zmax:
            zmax = z[k]
    if zmax >= h_i:
        return 0.0, 0.0, 0.0
    hc_in = (h_i - z[0]) / (1.5 + 0.5 * ke)
    q_in = w[0] * np.sqrt(G) * hc_in ** 1.5
    fb = _req(q_in, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out) - h_i
    if fb <= 0.0:
        n_l = n - 1
        q = q_in
        e_exit = z[n_l] + h_out[n_l] + q * q / (2.0 * G * (w[n_l] * h_out[n_l]) ** 2)
        diag[1] = abs(fb); diag[3] = 1.0
        return q, e_exit, fb + h_i
    # lower bracket end: R is NOT monotone as Q -> 0 at a free outfall (the brink critical depth becomes micrometric, the
    # Manning slope there is >> 1 and H_req jumps up). Halve Q from q_in until R < 0; flows below 2^-40 q_in are reported as 0.
    a = q_in
    fa = fb
    for _ in range(40):
        a = 0.5 * a
        fa = _req(a, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out) - h_i
        if fa < 0.0:
            break
    if fa >= 0.0:
        return 0.0, 0.0, 0.0
    b = q_in
    fb = _req(b, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out) - h_i
    c, fc = a, fa
    d = e = b - a
    tol_r = 1e-14 * (1.0 + h_i - z[0])
    it = 0
    for it in range(1, 100):
        if (fb > 0.0) == (fc > 0.0):
            c, fc = a, fa
            d = e = b - a
        if abs(fc) < abs(fb):
            a, b, c = b, c, b
            fa, fb, fc = fb, fc, fb
        tol = 2e-15 * abs(b) + 1e-300
        m = 0.5 * (c - b)
        if abs(m) <= tol or abs(fb) <= tol_r:
            break
        if abs(e) >= tol and abs(fa) > abs(fb):
            s = fb / fa
            if a == c:
                p = 2.0 * m * s
                qq = 1.0 - s
            else:
                qa = fa / fc
                r = fb / fc
                p = s * (2.0 * m * qa * (qa - r) - (b - a) * (r - 1.0))
                qq = (qa - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0.0:
                qq = -qq
            p = abs(p)
            if 2.0 * p < min(3.0 * m * qq - abs(tol * qq), abs(e * qq)):
                e = d
                d = p / qq
            else:
                d = m
                e = m
        else:
            d = m
            e = m
        a, fa = b, fb
        b = b + d if abs(d) > tol else b + (tol if m > 0 else -tol)
        fb = _req(b, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out) - h_i
    q = b
    # final profile at the returned discharge (h_out, hf_out, choked_out correspond to the last evaluation = b)
    req = fb + h_i
    n_l = n - 1
    e_exit = z[n_l] + h_out[n_l] + q * q / (2.0 * G * (w[n_l] * h_out[n_l]) ** 2)
    diag[0] = it; diag[1] = abs(fb); diag[2] = abs(c - b) / max(q, 1e-300); diag[3] = 2.0
    return q, e_exit, req


@njit(cache=True)
def edge_flow_diag(h_i, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out, diag):
    """Signed discharge (positive i -> j) and head bookkeeping; diag as in edge_flow_forward. Returns (Q, direction,
    E_entry, E_exit); h_out/hf_out/choked_out are indexed in the FLOW direction (index 0 = upstream end of the water)."""
    if h_i >= h_j:
        q, e_exit, e_entry = edge_flow_forward(h_i, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out, diag)
        return q, 1, e_entry, e_exit
    zr = z[::-1].copy()
    wr = w[::-1].copy()
    q, e_exit, e_entry = edge_flow_forward(h_j, h_i, zr, wr, n_man, dx, ke, msub, h_out, hf_out, choked_out, diag)
    return -q, -1, e_entry, e_exit


@njit(cache=True)
def edge_flow(h_i, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out):
    diag = np.zeros(4)
    return edge_flow_diag(h_i, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out, diag)


@njit(cache=True)
def smooth_activation(q, w0, h_floor):
    """Smooth switch-on of a discharge whose entrance critical depth h_c = (q^2/(g w0^2))^(1/3) is below h_floor.
    The raw Manning/critical-brink law has a fold (saddle-node) at micrometric depths: no root exists for heads below roughly 1e-5 m over
    the sill, then a root appears at Q ~ 1e-8 m3/s with h_c of a few micrometres, a regime in which Manning friction (a fully rough
    turbulent closure) is meaningless. A jump of that size is invisible to the water budget but makes an adaptive integrator
    stall at the moment a basin first overtops its sill. S(x) = x^2 (3 - 2x), x = h_c/h_floor, is 1 for h_c >= h_floor, so discharge
    with a critical depth of at least h_floor (1 mm by default) is untouched. h_floor = 0 returns the raw law."""
    if h_floor <= 0.0 or q <= 0.0:
        return q
    hc = (q * q / (G * w0 * w0)) ** (1.0 / 3.0)
    x = hc / h_floor
    if x >= 1.0:
        return q
    return q * x * x * (3.0 - 2.0 * x)
