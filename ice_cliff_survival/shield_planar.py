"""Resolved static two-material planar conduction vs the massless Robin shield (M2 Part 2).
Problem: debris layer -d<x<0 (k_d, rho c_d, alpha_d), initial T=0 (= T_m everywhere), hot-side face x=-d held at T_w from t=0 (liquid side assumed well mixed/fixed), ice boundary x=0 held at T_m=0 (a melting
interface stays at T_m); fixed geometry (interface motion and thickness change neglected). Heat delivered to the ice per unit area:
  resolved  q(t) = (k_d T_w/d) [1 + 2 sum_{n>=1} (-1)^n exp(-n^2 pi^2 alpha t/d^2)],   Q(t) = int q
  massless  q = k_d T_w/d   (the M1 shield h=k/d), Q_m = k_d T_w t/d
  deficit   Q_m - Q -> rho c d T_w/6 = t_lag * q_inf   with  t_lag = d^2/(6 alpha)  [exact late-time limit of the series: sum (-1)^n/n^2 = -pi^2/12]
The resolved reference is checked against an independent implicit finite-volume solution (fine grid, Crank-Nicolson) in the tests."""
import math
import numpy as np
def q_series(t, d, k, rhoc, Tw, nterms=400):
    a = k / rhoc; n = np.arange(1, nterms + 1); t = np.asarray(t, float)[..., None]
    return k * Tw / d * (1.0 + 2.0 * np.sum((-1.0) ** n * np.exp(-(n * math.pi) ** 2 * a * t / d ** 2), axis=-1))
def Q_series(t, d, k, rhoc, Tw, nterms=400):
    a = k / rhoc; n = np.arange(1, nterms + 1); t = np.asarray(t, float)[..., None]; lam = (n * math.pi) ** 2 * a / d ** 2
    return k * Tw / d * (t[..., 0] + 2.0 * np.sum((-1.0) ** n * (1 - np.exp(-lam * t)) / lam, axis=-1))
def Q_massless(t, d, k, Tw): return k * Tw / d * np.asarray(t, float)
def t_lag(d, k, rhoc): return d * d * rhoc / (6.0 * k)
def fv_numeric(t_end, d, k, rhoc, Tw, N=800, nsteps=4000):
    """independent: cell-centred FV, Crank-Nicolson after 20 backward-Euler start-up steps (Rannacher), Dirichlet T=Tw at x=-d, T=0 at x=0; returns cumulative heat into the ice boundary (flux x energy)."""
    from scipy.linalg import solve_banded
    dx = d / N; a = k / rhoc; T = np.zeros(N)       # x from hot (index 0) to ice (N-1)
    dt = t_end / nsteps; Q = 0.0; r = a * dt / dx ** 2
    def Lop(T):
        Tg = np.concatenate(([2 * Tw - T[0]], T, [-T[-1]])); return (Tg[:-2] - 2 * Tg[1:-1] + Tg[2:]) / dx ** 2 * a      # ghost cells give Dirichlet faces
    def mats(theta):
        ab = np.zeros((3, N)); ab[1] = 1 + theta * a * dt / dx ** 2 * 2; ab[0, 1:] = -theta * a * dt / dx ** 2; ab[2, :-1] = -theta * a * dt / dx ** 2; ab[1, 0] += theta * a * dt / dx ** 2; ab[1, -1] += theta * a * dt / dx ** 2; return ab
    def flux_ice(T): return k * (0 - T[-1]) / (dx / 2) * -1      # heat into the ice boundary = -k dT/dx toward ice -> k T[-1]/(dx/2)
    for n in range(nsteps):
        th = 1.0 if n < 20 else 0.5; ab = mats(th)
        # Dirichlet source enters through the ghost cell: rhs = T + (1-th) dt L(T) + th dt a/dx^2 2 Tw e_0
        rhs = T + (1 - th) * dt * Lop(T) + th * dt * a / dx ** 2 * 2 * Tw * np.eye(1, N, 0)[0]
        Tn = solve_banded((1, 1), ab, rhs); Q += 0.5 * dt * (k * T[-1] / (dx / 2) + k * Tn[-1] / (dx / 2)) if th == 0.5 else dt * k * Tn[-1] / (dx / 2); T = Tn
    return Q

def study(out):
    import json
    rows = []
    for d in (0.002, 0.007, 0.02, 0.1):
        for rhoc in (1.5e6, 2.5e6):
            for k in (0.3, 0.5, 1.0):
                Tw = 0.00631 * 3.34e5 / 4217.0; tl = t_lag(d, k, rhoc)
                for t in (60.0, 600.0, 3600.0, 86400.0):
                    Q, Qm = float(Q_series(t, d, k, rhoc, Tw)), float(Q_massless(t, d, k, Tw)); rows.append(dict(d=d, rhoc=rhoc, k=k, t=t, t_lag=tl, rel_deficit=(Qm - Q) / Qm, rel_deficit_asymptotic=tl / t))
    json.dump(dict(rows=rows, note="resolved = analytic series (checked against independent FV in tests); fixed geometry, T_w step, ice face at T_m; ASSUMED inputs, not measurements"), open(out, "w"), indent=1); return rows
if __name__ == "__main__":
    import sys; study(sys.argv[1])
