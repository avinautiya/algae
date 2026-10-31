"""Separate one-dimensional Stefan verification (NOT OpenFOAM): one-phase melting, liquid on 0<x<s(t) held at T_w at x=0, ice at the melting point beyond s (so no conduction in ice).
Exact Neumann solution: s = 2 lam sqrt(kappa t), lam exp(lam^2) erf(lam) = Ste/sqrt(pi), Ste = c (T_w - T_m)/L.   Numerics: enthalpy method on a fixed grid (apparent heat capacity-free, exact enthalpy conservation),
three resolutions and a separate time-step refinement; front position from the liquid fraction. Equal densities (the density jump enters only through the normal-velocity condition u_w.n = (rho_i/rho_w - 1) m_n, see MATHEMATICAL_UPGRADES.md section 3, not tested here)."""
import json, os
import numpy as np
from scipy.special import erf
from scipy.optimize import brentq
Ste = 0.2; kappa = 1.0; Lc = 1.0; tend = 1.0
lam = brentq(lambda l: l * np.exp(l * l) * erf(l) - Ste / np.sqrt(np.pi), 1e-6, 5.0); s_exact = 2 * lam * np.sqrt(kappa * tend)


def run(nx, dt):
    """cells on [0, X]; unknown: enthalpy per cell e = T (liquid) or latent fraction f in [0,1] (at T=0 partially melted), nondimensional c=1, L=1/Ste, T_w = 1."""
    X = 1.0; dx = X / nx; Lh = 1.0 / Ste; Hn = np.zeros(nx)          # H = c T + L f ; H in [0, L] mixed (T=0), H>L liquid, T = H - L
    t = 0.0; n = int(round(tend / dt))
    for _ in range(n):
        T = np.where(Hn > Lh, Hn - Lh, 0.0); Tw = np.concatenate([[1.0], T, [0.0]])        # wall T_w=1; far end at T_m
        # wall cell: half-cell conduction from the wall value, conservative explicit-in-time update via implicit tridiagonal for liquid is avoided: use small dt (explicit, stable dt <= dx^2/2)
        flux = -(Tw[1:] - Tw[:-1]) / dx; flux[0] = -(T[0] - 1.0) / (0.5 * dx)
        Hn = Hn - dt / dx * (flux[1:] - flux[:-1]); t += dt
    f = np.clip(Hn / Lh, 0.0, 1.0); s = float(np.sum(f) * dx)          # melted length = integral of liquid fraction, exact for the sharp front
    return s, float(np.sum(Hn) * dx)


out = []
for nx in (50, 100, 200):
    dt = 0.2 * (1.0 / nx) ** 2; s, E = run(nx, dt); out.append(dict(nx=nx, dt=dt, s=s, err=abs(s - s_exact), rel=abs(s - s_exact) / s_exact))
for a, b in zip(out[:-1], out[1:]): b["order"] = float(np.log2(a["err"] / b["err"]))
tt = []
for dtf in (1.0, 0.5, 0.25):      # separate time-step refinement at nx=100 (explicit limit dt = 0.2 dx^2 scaled)
    dt = 0.2 * (1.0 / 100) ** 2 * dtf; s, E = run(100, dt); tt.append(dict(dt=dt, s=s, err=abs(s - s_exact)))
res = dict(Ste=Ste, lam=lam, s_exact=s_exact, space=out, time=tt, note="explicit enthalpy scheme: time error small relative to space error; order ~1 expected for front position")
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "stefan1d_results.json"), "w"), indent=1); print(json.dumps(res, indent=1))
