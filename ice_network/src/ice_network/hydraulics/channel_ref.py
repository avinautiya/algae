"""REFERENCE (slow, bisection-only) implementation of the discharge calculation. Kept unchanged as the baseline against which the fast
solver in channel.py is verified; do not use in production.

Quasi-steady open-channel hydraulics of one reach connecting two storage nodes.

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
def _crit_depth_ref(q, w):
    return (q * q / (G * w * w)) ** (1.0 / 3.0)


@njit(cache=True)
def _sf_ref(q, w, h, n_man):
    a = w * h
    r = a / (w + 2.0 * h)
    return (n_man * q / (a * r ** (2.0 / 3.0))) ** 2


@njit(cache=True)
def _station_depth_ref(e_next, sf_next, z, w, q, n_man, dx):
    """Implicit standard-step station: find the subcritical depth h with
        z + h + q^2/(2 g w^2 h^2) = e_next + 0.5 (Sf(h) + sf_next) dx
    by bisection on [h_c, h_hi] (robust: no fixed-point iteration, which oscillates near critical flow).
    Returns (h, E_k, choked). If no subcritical root exists the depth is critical and the energy line is raised
    to z + 1.5 h_c (control point)."""
    q2 = q * q / (w * w)
    hc = (q2 / G) ** (1.0 / 3.0)
    sf_c = _sf_ref(q, w, hc, n_man)
    g_c = z + hc + q2 / (2.0 * G * hc * hc) - e_next - 0.5 * (sf_c + sf_next) * dx
    if g_c >= 0.0:
        return hc, z + 1.5 * hc, True
    lo = hc
    hi = e_next + 0.5 * (sf_c + sf_next) * dx - z
    if hi <= lo:
        hi = 2.0 * lo
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        gm = z + mid + q2 / (2.0 * G * mid * mid) - e_next - 0.5 * (_sf_ref(q, w, mid, n_man) + sf_next) * dx
        if gm > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-15 * hi:
            break
    h = 0.5 * (lo + hi)
    return h, e_next + 0.5 * (_sf_ref(q, w, h, n_man) + sf_next) * dx, False


@njit(cache=True)
def gvf_upstream_ref(q, h_l, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out):
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
    sf_next = _sf_ref(q, w[n - 1], h_l, n_man)
    dxs = dx / msub
    h_k = h_l
    for k in range(n - 2, -1, -1):
        hf = 0.0
        ch_any = 0
        for j in range(1, msub + 1):
            f = j / msub
            zz = z[k + 1] + (z[k] - z[k + 1]) * f
            ww = w[k + 1] + (w[k] - w[k + 1]) * f
            h_k, e_k, ch = _station_depth_ref(e_next, sf_next, zz, ww, q, n_man, dxs)
            sf_k = _sf_ref(q, ww, h_k, n_man)
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
def edge_flow_forward_ref(h_i, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out):
    """Discharge from node i to node j for h_i >= h_j. Returns (Q, E_exit_total, E_entry_total). Q >= 0."""
    n = z.size
    if h_i <= z[0] or h_i <= h_j:
        return 0.0, 0.0, 0.0
    # inlet-control capacity
    hc_in = (h_i - z[0]) / (1.5 + 0.5 * ke)
    q_in = w[0] * np.sqrt(G) * hc_in ** 1.5
    # dry-limit check: a hump higher than the reservoir stage blocks the channel
    zmax = z[0]
    for k in range(n):
        if z[k] > zmax:
            zmax = z[k]
    if zmax >= h_i:
        return 0.0, 0.0, 0.0
    q_lo, q_hi = 0.0, q_in

    def_req = 0.0
    # R(Q) = H_req(Q) - h_i, increasing in Q
    # evaluate at q_hi
    hl = max(h_j - z[n - 1], _crit_depth_ref(q_hi, w[n - 1]))
    req = gvf_upstream_ref(q_hi, hl, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out)
    if req <= h_i:
        q = q_hi
        e_exit = z[n - 1] + h_out[n - 1] + q * q / (2.0 * G * (w[n - 1] * h_out[n - 1]) ** 2)
        return q, e_exit, req
    for _ in range(70):
        q_mid = 0.5 * (q_lo + q_hi)
        if q_mid <= 0.0:
            break
        hl = max(h_j - z[n - 1], _crit_depth_ref(q_mid, w[n - 1]))
        req = gvf_upstream_ref(q_mid, hl, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out)
        if req > h_i:
            q_hi = q_mid
        else:
            q_lo = q_mid
        if q_hi - q_lo < 1e-14 * max(q_hi, 1e-12):
            break
    q = q_lo
    if q <= 0.0:
        return 0.0, 0.0, 0.0
    hl = max(h_j - z[n - 1], _crit_depth_ref(q, w[n - 1]))
    req = gvf_upstream_ref(q, hl, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out)
    e_exit = z[n - 1] + h_out[n - 1] + q * q / (2.0 * G * (w[n - 1] * h_out[n - 1]) ** 2)
    return q, e_exit, req


@njit(cache=True)
def edge_flow_ref(h_i, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out):
    """Signed discharge (positive i -> j) and head bookkeeping. Returns (Q, direction, E_entry, E_exit) where the head
    quantities refer to the direction of flow; h_out/hf_out/choked_out are indexed in the FLOW direction
    (index 0 = upstream end of the flowing water)."""
    if h_i >= h_j:
        q, e_exit, e_entry = edge_flow_forward_ref(h_i, h_j, z, w, n_man, dx, ke, msub, h_out, hf_out, choked_out)
        return q, 1, e_entry, e_exit
    zr = z[::-1].copy()
    wr = w[::-1].copy()
    q, e_exit, e_entry = edge_flow_forward_ref(h_j, h_i, zr, wr, n_man, dx, ke, msub, h_out, hf_out, choked_out)
    return -q, -1, e_entry, e_exit
