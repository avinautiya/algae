"""Water energy balance along a reach and heat delivered to the ice wall (temperate ice, wall at T_m = 0 C).

Per segment k (between stations k and k+1; mean width w, mean depth h, wetted perimeter P = w + 2h):
    rho_w c_w Q dT/dx = phi' - hT P (T - T_m)           phi' = dissipation heat per length (uniform in the segment)
with h_T = Nu k_w / D_h, Nu = C_nu Re^0.8 Pr^0.3 (Dittus-Boelter, water being cooled), D_h = 4 A / P.
Over a segment with constant coefficients the solution is exact:
    T_{k+1} = T_eq + (T_k - T_eq) exp(-hT P dx / (rho_w c_w Q)),   T_eq = T_m + phi' / (hT P).
Heat delivered to the ice in the segment (W): q = rho_w c_w Q (T_k - T_{k+1}) + Phi_seg  (>= 0 for T >= T_m),
so the segment energy balance closes by construction; ice melted: dV_ice/dt = q / (rho_i L_f).
Dissipation is an EXPLICIT heat source here; the melt law contains no separate incision term (no double counting).
"""
import numpy as np
from numba import njit

from ..constants import CP_W, K_W, NU_W, PR_W, RHO_W


@njit(cache=True)
def heat_transfer_coeff(q, w, h, nu_factor):
    a = w * h
    p = w + 2.0 * h
    dh = 4.0 * a / p
    v = q / a
    re = max(v * dh / NU_W, 1.0)
    nu = 0.023 * re ** 0.8 * PR_W ** 0.3
    return nu_factor * nu * K_W / dh


@njit(cache=True)
def reach_heat(q, t_in, t_m, h, w, dx, phi_seg, phi_inlet, nu_factor, t_out_arr, q_ice):
    """March the water temperature down the reach (arrays in flow direction).
    phi_seg[k]: friction dissipation in segment k (W); phi_inlet: entrance-loss dissipation (W), added to segment 0.
    Fills t_out_arr[k] (temperature at station k, length N) and q_ice[k] (W delivered to ice in segment k).
    Returns the temperature at the downstream end."""
    n = h.size
    t_out_arr[0] = t_in
    rc = RHO_W * CP_W
    for k in range(n - 1):
        wm = 0.5 * (w[k] + w[k + 1])
        hm = 0.5 * (h[k] + h[k + 1])
        p = wm + 2.0 * hm
        htc = heat_transfer_coeff(q, wm, hm, nu_factor)
        phi = phi_seg[k] + (phi_inlet if k == 0 else 0.0)
        # equilibrium excess temperature and exact exponential relaxation
        hp = htc * p
        t_eq = t_m + (phi / dx) / hp
        decay = np.exp(-hp * dx / (rc * q))
        t_next = t_eq + (t_out_arr[k] - t_eq) * decay
        t_out_arr[k + 1] = t_next
        q_ice[k] = rc * q * (t_out_arr[k] - t_next) + phi
    return t_out_arr[n - 1]
