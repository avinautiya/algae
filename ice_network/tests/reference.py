"""Independent reference solutions used only by the tests (different formulation from the solver).

gvf_ode_Q: steady gradually varied flow in a PRISMATIC rectangular channel solved as the depth ODE
    dh/dx = (S0 - Sf) / (1 - Fr^2)       (x measured downstream, S0 = bed slope, Sf Manning friction slope)
integrated upstream with a high-accuracy adaptive Runge-Kutta, and the discharge found by matching the reservoir
stage at the entrance, H_i = z_0 + h_0 + (1 + k_e) v_0^2 / (2 g). The solver under test marches the TOTAL HEAD
equation on a station grid, so agreement is a check of the numerics and of the implementation, not of the closure
(Manning, entrance loss), which both share.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G = 9.81


def sf_manning(q_tot, w, h, n):
    a = w * h
    r = a / (w + 2 * h)
    return (n * q_tot / (a * r ** (2 / 3))) ** 2


def profile_ode(Q, w, z0, z1, L, n, h_down, rtol=1e-11):
    s0 = (z0 - z1) / L
    hc = (Q * Q / (G * w * w)) ** (1 / 3)

    def f(x, h):
        hh = h[0]
        fr2 = Q * Q / (G * w * w * hh ** 3)
        return [(s0 - sf_manning(Q, w, hh, n)) / (1 - fr2)]

    # integrate from x = L (downstream) to x = 0 (upstream): use s = L - x
    sol = solve_ivp(lambda s, h: [-f(L - s, h)[0]], (0, L), [h_down], method="DOP853", rtol=rtol, atol=1e-13, dense_output=True)
    return sol


def gvf_ode_Q(h_res, h_tail, w, z0, z1, L, n, ke):
    """Discharge from a reservoir at stage h_res through a prismatic rectangular channel to a tailwater stage."""
    def required(Q):
        hc = (Q * Q / (G * w * w)) ** (1 / 3)
        h_down = max(h_tail - z1, hc * (1 + 1e-7))
        sol = profile_ode(Q, w, z0, z1, L, n, h_down)
        h0 = sol.y[0, -1]
        v0 = Q / (w * h0)
        return z0 + h0 + (1 + ke) * v0 * v0 / (2 * G)
    q_in = w * np.sqrt(G) * ((h_res - z0) / (1.5 + 0.5 * ke)) ** 1.5
    return brentq(lambda Q: required(Q) - h_res, 1e-9, q_in * (1 - 1e-9), xtol=1e-14, rtol=1e-13)
