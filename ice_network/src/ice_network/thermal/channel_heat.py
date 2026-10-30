"""Water energy balance along a reach and heat delivered to the ice wall (temperate ice, wall at T_m = 0 C).

Per segment k (between stations k and k+1; mean width w, mean depth h, wetted perimeter P = w + 2h):
    rho_w c_w Q dT/dx = phi' - hT P (T - T_m)           phi' = dissipation heat per length (uniform in the segment)
with h_T = Nu k_w / D_h, Nu = C_nu Re^0.8 Pr^0.3 (Dittus-Boelter, water being cooled), D_h = 4 A / P.
Over a segment with constant coefficients the solution is exact. It is evaluated in the cancellation-free form of
`segment_heat` (a = h_T P dx / (rho_w c_w Q), F = 1 - exp(-a)); the earlier equilibrium-temperature form
T_eq = T_m + phi'/(h_T P) lost all digits when a ~ 1e-13 (stress test nu_factor = 1e-12; see
ice_network_verification/records/claude_heat_kernel_audit.json). That was a numerical stress condition, not a default setting.
Heat delivered to the ice in the segment (W): q = rho_w c_w Q (T_k - T_{k+1}) + Phi_seg, so the segment energy balance closes
by construction; ice melted: dV_ice/dt = q / (rho_i L_f).
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
def ice_share_of_source(a):
    """1 - (1 - exp(-a))/a = (a - 1 + exp(-a))/a, evaluated without cancellation for small a (series for a < 0.5)."""
    if a < 0.5:
        # sum_{k>=1} (-1)^(k+1) a^k / (k+1)!   (a/2 - a^2/6 + a^3/24 - ...)
        term = 0.5 * a
        total = term
        for k in range(2, 16):
            term *= -a / (k + 1)
            total += term
        return total
    return 1.0 - (-np.expm1(-a)) / a


@njit(cache=True)
def segment_heat(c_rate, ua, theta_in, phi):
    """Exact solution of C dT/dx-type balance over one segment with constant coefficients (stable algebra).
    C = rho_w c_w Q (W/K), UA = h_T P dx (W/K), theta_in = T_in - T_m (K), phi = dissipation heated into the water in
    the segment (W, uniform). With a = UA/C and F = 1 - exp(-a):
        theta_out = theta_in exp(-a) + (phi/C) (F/a)
        q_ice     = C theta_in F + phi (1 - F/a)
    F uses expm1; (1 - F/a) uses the series for small a; at a = 0 theta_out = theta_in + phi/C and q_ice = 0.
    No clipping: a negative theta_in produces a negative q_ice (refreezing), which this melt-only model does not support.
    Returns (theta_out, q_ice)."""
    a = ua / c_rate
    f = -np.expm1(-a)
    if a < 1e-300:
        return theta_in + phi / c_rate, 0.0
    share_ice = ice_share_of_source(a)
    theta_out = theta_in * np.exp(-a) + (phi / c_rate) * (1.0 - share_ice)
    q_ice = c_rate * theta_in * f + phi * share_ice
    return theta_out, q_ice


@njit(cache=True)
def reach_heat(q, t_in, t_m, h, w, dx, phi_seg, phi_inlet, nu_factor, t_out_arr, q_ice):
    """March the water temperature down the reach (arrays in flow direction), segment by segment (segment lengths
    may differ only through dx; geometry may vary).
    phi_seg[k]: friction dissipation in segment k (W); phi_inlet: entrance-loss dissipation (W), added to segment 0.
    Fills t_out_arr[k] (temperature at station k, length N) and q_ice[k] (W delivered to ice in segment k).
    Requires q > 0. Returns the temperature at the downstream end."""
    n = h.size
    t_out_arr[0] = t_in
    rc = RHO_W * CP_W
    c_rate = rc * q
    for k in range(n - 1):
        wm = 0.5 * (w[k] + w[k + 1])
        hm = 0.5 * (h[k] + h[k + 1])
        p = wm + 2.0 * hm
        htc = heat_transfer_coeff(q, wm, hm, nu_factor)
        phi = phi_seg[k] + (phi_inlet if k == 0 else 0.0)
        th_out, qk = segment_heat(c_rate, htc * p * dx, t_out_arr[k] - t_m, phi)
        t_out_arr[k + 1] = t_m + th_out
        q_ice[k] = qk
    return t_out_arr[n - 1]
