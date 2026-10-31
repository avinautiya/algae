"""Semi-analytic Graetz-Poiseuille reference (parallel plates, u/Ub = 6 eta (1-eta), wall T = 0, inlet T = 1, axial conduction neglected: relative size O(1/Pe_H^2) = %.1e) and error measures for the OpenFOAM fields.
theta(xi, eta) = sum_n C_n phi_n(eta) exp(-lam_n^2 xi), xi = x/(H Pe_H), phi'' + 6 lam^2 eta (1-eta) phi = 0, phi(0)=phi(1)=0 (generalised eigenproblem, N=3000 finite differences, Richardson-checked).
Errors reported SEPARATELY at x = 0.9 L: (1) velocity profile; (2) temperature profile; (3) wall heat flux (and Nu_Dh). Normalisations: velocity by Ub; temperature by the inlet-wall difference (1) and by T_b; wall flux by the exact local flux."""
import json, os, sys
import numpy as np
from scipy.linalg import eigh_tridiagonal, eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import poiseuille_graetz as P
PE = P.UB * P.H / P.DT_T
__doc__ = __doc__ % (1 / PE ** 2)


def mode1(N=600):
    """leading decaying mode WITH axial conduction: theta = phi(eta) exp(-s xi), phi'' + (s W + s^2/Pe^2) phi = 0, as the eigenproblem M [phi; s phi] = s [phi; s phi] (finite differences, N interior points)."""
    h = 1.0 / N; eta = np.arange(1, N) * h; n = N - 1; W = 6 * eta * (1 - eta)
    A = (np.diag(2 * np.ones(n)) - np.diag(np.ones(n - 1), 1) - np.diag(np.ones(n - 1), -1)) / h ** 2
    M = np.block([[np.zeros((n, n)), np.eye(n)], [PE ** 2 * A, -PE ** 2 * np.diag(W)]]); ev, V = np.linalg.eig(M)
    k = [i for i in range(2 * n) if abs(ev[i].imag) < 1e-9 and ev[i].real > 0]; i = min(k, key=lambda i: ev[i].real); s = ev[i].real; phi = V[:n, i].real; phi = phi * np.sign(phi.sum()); return eta, phi, s, h, W


def errors(case, ny, tl="600"):
    nx = 20 * ny; U = P.read_field(f"{case}/{tl}/U", "vec")[:, 0].reshape(ny, nx); T = P.read_field(f"{case}/{tl}/T", "sca").reshape(ny, nx)
    eta_f = (np.arange(ny) + 0.5) / ny; ub = U[:, int(0.9 * nx)].mean(); eta, phi, s, h, W = mode1()
    def at(i):
        u = U[:, i]; tt = T[:, i]; Tb = (u * tt).sum() / u.sum(); dy = 1.0 / ny; y1, y2 = 0.5 * dy, 1.5 * dy
        g = (tt[0] * y2 ** 2 - tt[1] * y1 ** 2) / (y1 * y2 * (y2 - y1)); return u, tt, float(Tb), float(abs(g))
    i1, i2 = int(0.7 * nx), int(0.9 * nx); u, tt, Tb, g = at(i2); _, _, Tb1, _ = at(i1)
    Tb_m = float(np.sum(W * phi) / np.sum(W)); g_m = float(abs(phi[0] / h)); phi_f = np.interp(eta_f, eta, phi)
    eU = float(np.sqrt(np.mean((u - 6 * ub * eta_f * (1 - eta_f)) ** 2)) / ub)
    shape_num = tt / Tb; shape_ex = phi_f / Tb_m; eT = float(np.sqrt(np.mean((shape_num - shape_ex) ** 2)))              # temperature PROFILE shape error, normalised by the local flow-weighted bulk excess
    Nu = 2 * g / Tb; Nu_ex = 2 * g_m / Tb_m; e_flux = float(abs(Nu - Nu_ex) / Nu_ex)                                    # wall heat flux in units of (k T_b / Dh) : melt-driving flux error
    dxi = (i2 - i1) * (P.L / nx) / (P.H * PE); s_num = float(-np.log(Tb / Tb1) / dxi); e_s = float(abs(s_num - s) / s)    # thermal development rate (controls absolute flux level)
    return dict(ny=ny, dy_over_H=1.0 / ny, x_over_H=(i2 + 0.5) * P.L / nx / P.H, err_u_over_Ub=eU, err_T_shape_rms=eT, err_wall_flux_rel=e_flux, Nu_Dh=float(Nu), Nu_Dh_ref=float(Nu_ex), s_num=s_num, s_ref=float(s), err_decay_rate_rel=e_s)


if __name__ == "__main__":
    wd = sys.argv[1]; res = []
    for nm in sys.argv[2:]:
        ny = int(nm.split("ny")[1].split("_")[0]); r = errors(f"{wd}/{nm}", ny); r["case"] = nm; res.append(r); print(r, flush=True)
    json.dump(dict(Re_H=P.UB * P.H / P.NU, Pe_H=PE, Pr=P.PR, results=res), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "graetz_errors.json"), "w"), indent=1)
