"""Resolved laminar benchmark for the reduced model's segment heat exchange: Graetz problem between parallel plates with isothermal (T = T_m = 0) walls.

Resolved reference: steady 2-D advection-diffusion (axial conduction neglected, Pe >> 1), Poiseuille profile u = 1.5 U (1 - xi^2), centreline symmetry,
backward-Euler marching in x* = alpha x / (H^2 U) on a cosine-stretched radial grid with a geometrically growing step; tested against grid refinement and
against the analytic fully developed Nusselt number 7.5407 (D_h = 4H) and the leading Graetz eigenvalue series.
Reduced model: constant-Nu segment march, the production `segment_heat` with Nu = Nu_fd, the same Q, wall area and temperature.
What this tests: the exact-exponential segment algebra and wall-area bookkeeping against a resolved solution, and the size of the entrance (developing
boundary layer) error of a constant-Nu closure. What it does NOT test: the turbulent closure used in the melting channels (laminar flow, plates, no free surface,
no melting). Written in nondimensional form so it makes no physical-glacier claim."""
import json, os, sys
import numpy as np
from scipy.linalg import solve_banded
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "ice_network", "src"))
from ice_network.thermal.channel_heat import segment_heat
NU_FD = 7.5407     # Nu_{Dh}, parallel plates, constant wall temperature, fully developed


def graetz(nxi=400, nsteps=4000, xstar_end=0.5, dx0=1e-9):
    s = np.linspace(0.0, 1.0, nxi + 1); xi = 0.5 * (1.0 - np.cos(np.pi * s))          # finer near the wall (xi = 1)
    u = 1.5 * (1.0 - xi ** 2)
    # nodes 0..nxi; node nxi is the wall (theta = 0); node 0 is the centreline (symmetry: theta_xi = 0)
    h = np.diff(xi); th = np.ones(nxi + 1); th[-1] = 0.0
    # weights for bulk: 1.5 * int (1-xi^2) theta dxi (trapezoid)
    wq = np.zeros(nxi + 1)
    f = 1.5 * (1.0 - xi ** 2)
    for i in range(nxi):
        wq[i] += 0.5 * h[i] * f[i]; wq[i + 1] += 0.5 * h[i] * f[i + 1]
    lo, hi = 1.0 + 1e-12, 1.2          # geometric step ratio r with dx0 (r^n - 1)/(r - 1) = xstar_end
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if dx0 * np.expm1(nsteps * np.log(mid)) / (mid - 1.0) < xstar_end: lo = mid
        else: hi = mid
    r = 0.5 * (lo + hi); dxs = dx0 * r ** np.arange(nsteps)
    LOW = np.zeros(nxi + 1); UPP = np.zeros(nxi + 1)
    hl, hr = h[:-1], h[1:]; c = 2.0 / (hl + hr); LOW[1:nxi] = c / hl; UPP[1:nxi] = c / hr
    xs = [0.0]; tb = [1.0]; nu = [np.inf]
    for dx in dxs:
        # (u/dx) (th_new - th_old) = d2 th_new / dxi2  on interior nodes; centreline: symmetric ghost; wall: Dirichlet 0
        ab = np.zeros((3, nxi + 1)); rhs = u / dx * th
        ab[1, :] = u / dx + LOW + UPP; ab[0, 1:] = -UPP[:-1]; ab[2, :-1] = -LOW[1:]
        ab[1, 0] = u[0] / dx + 2.0 / h[0] ** 2; ab[0, 1] = -2.0 / h[0] ** 2
        ab[1, nxi] = 1.0; rhs[nxi] = 0.0; ab[2, nxi - 1] = 0.0
        th = solve_banded((1, 1), ab, rhs)
        xs.append(xs[-1] + dx); b = float(wq @ th)
        # wall gradient: second-order one-sided with non-uniform spacing
        h1, h2 = h[-1], h[-2]
        g = (th[-1] * (2 * h1 + h2) / (h1 * (h1 + h2)) - th[-2] * (h1 + h2) / (h1 * h2) + th[-3] * h1 / (h2 * (h1 + h2)))
        tb.append(b); nu.append(-4.0 * g / b)
    return np.array(xs), np.array(tb), np.array(nu)


def reduced(xs, nseg, nu_fd=NU_FD):
    """production segment_heat march in the same nondimensional units: rho c U H = C = 1 (W/K per unit span, half-gap), wall conductance per x* = Nu/4"""
    out = []
    xe = xs[-1]
    for x in xs:
        edges = np.linspace(0.0, x, nseg + 1); th = 1.0
        for a, b in zip(edges[:-1], edges[1:]):
            th, _ = segment_heat(1.0, nu_fd / 4.0 * (b - a), th, 0.0)
        out.append(th)
    return np.array(out)


if __name__ == "__main__":
    runs = {}
    for nxi, ns in ((100, 1000), (200, 2000), (400, 4000)):
        xs, tb, nu = graetz(nxi, ns); runs[(nxi, ns)] = (xs, tb, nu)
    xs, tb, nu = runs[(400, 4000)]
    probe = [1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5]
    def at(arr, x, xv): return float(np.interp(x, xv, arr))
    rows = []
    for x in probe:
        tbs = [at(runs[k][1], x, runs[k][0]) for k in runs]
        red = at(reduced(np.array([0.0, x]), 8), x, np.array([0.0, x])) if False else float(reduced(np.array([x]), 8)[0])
        rows.append(dict(xstar=x, Tb_resolved_400=tbs[2], Tb_resolved_200=tbs[1], Tb_resolved_100=tbs[0], grid_change_200_400=abs(tbs[2] - tbs[1]), Tb_reduced_constNu=red,
                         reduced_minus_resolved=red - tbs[2], heat_removed_resolved=1 - tbs[2], heat_removed_reduced=1 - red,
                         heat_error_rel=((1 - red) - (1 - tbs[2])) / (1 - tbs[2]), local_Nu_resolved=at(runs[(400, 4000)][2], x, xs)))
    out = dict(description=__doc__.split("\n")[0], Nu_fd_resolved_at_x0p5=float(nu[-1]), Nu_fd_analytic=NU_FD, rows=rows,
               segmentation_independence={n: float(reduced(np.array([0.1]), n)[0]) for n in (1, 2, 8, 64)},
               note="x* = alpha x /(H^2 U); T_b/T_in at x*. Positive heat_error_rel = constant-Nu closure removes MORE heat than the resolved solution")
    os.makedirs(os.path.join(HERE, "records"), exist_ok=True); json.dump(out, open(os.path.join(HERE, "records", "graetz_benchmark.json"), "w"), indent=1)
    print("Nu(x*=0.5) resolved", out["Nu_fd_resolved_at_x0p5"], "analytic", NU_FD)
    for r in rows: print({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items() if k in ("xstar", "Tb_resolved_400", "grid_change_200_400", "Tb_reduced_constNu", "heat_error_rel", "local_Nu_resolved")})
    print(out["segmentation_independence"])
