"""Classical and local baselines for the matched-state experiment (frozen in FREEZE_TH.md).
Film profile w(eta) = 1.5 (2 eta - eta^2) (mean 1; no-slip cold bottom at eta = 0, zero-shear adiabatic top at eta = 1); theta_xi w = theta_eta_eta, xi = (x - x0)/(H Pe_H).
Graetz series (axial conduction neglected; Jones 1975 shows it matters only within O(1/Pe) of the wall-temperature step; here Pe_H = 350)  theta(xi, eta) = sum C_n phi_n e^{-lam_n^2 xi}, phi'' + lam^2 w phi = 0, phi(0) = 0, phi'(1) = 0.
Nu_H := (1/H) q_w H/(k T_b) defined by q_w = (k/H) Nu_H T_b with T_b = int w theta; Nu_fd from the first mode."""
import numpy as np
from scipy.linalg import eigh
PE = 5e-3 * 0.01 / (1e-6 / 7.0)
def modes(N=2000, nm=120):
    h = 1.0 / N; eta = (np.arange(N) + 0.5) * h; w = 1.5 * (2 * eta - eta ** 2)
    A = (np.diag(2 * np.ones(N)) - np.diag(np.ones(N - 1), 1) - np.diag(np.ones(N - 1), -1)); A[0, 0] = 3.0; A[-1, -1] = 1.0; A = A / h ** 2
    lam2, V = eigh(A, np.diag(w), subset_by_index=[0, nm - 1]); phi = V / np.sqrt(np.sum(w[:, None] * V ** 2, axis=0) * h)[None, :]    # w-normalised
    C = np.array([np.sum(w * phi[:, k]) * h for k in range(nm)])          # unit inlet temperature: C_n = int w phi_n (phi_n w-orthonormal)
    dphi0 = phi[0, :] * 2.0 / h        # derivative at wall (phi_0 = -phi_1 ghost => slope = 2 phi_1/h)... cell-centred: phi(h/2)=phi_1, phi(0)=0 -> first-order; refined below
    # second-order wall slope from the two nearest cells (quadratic through 0, phi_1, phi_2)
    x1, x2 = 0.5 * h, 1.5 * h; dphi0 = (phi[0, :] * x2 ** 2 - phi[1, :] * x1 ** 2) / (x1 * x2 * (x2 - x1))
    return eta, w, lam2, phi, C, dphi0, h
_M = None
def M():
    global _M
    if _M is None: _M = modes()
    return _M
def graetz(xi):
    eta, w, lam2, phi, C, dphi0, h = M(); e = np.exp(-lam2 * xi) * C; theta = phi @ e; Tb = float(np.sum(w * theta) * h); qw = float(np.sum(dphi0 * e))      # d theta/d eta at wall
    return eta, theta, Tb, qw
def Nu_fd():
    eta, w, lam2, phi, C, dphi0, h = M(); Tb1 = float(np.sum(w * phi[:, 0]) * h); return float(dphi0[0] / Tb1), float(lam2[0])
def graetz_section(x_from_start_H):
    """bulk T, wall flux q = (alpha/H) dtheta/deta (units of alpha/H, unit inlet) at a distance x (in H) downstream of the cold-wall start"""
    xi = x_from_start_H / PE; eta, theta, Tb, qw = graetz(xi); return Tb, qw, theta, eta
def local_closure(Tb_c, kind="step", X=30.0, n=3000):
    """history-blind marching: dTb/dx = -(alpha Nu/(U_loc H_loc^2)) Tb ; returns x [H], q [alpha/H units, per H_loc basis -> physical q/(alpha/H)] and Q(X) integral (units alpha*H*..., same as the CFD observable per unit depth with lengths in H)."""
    nu, _ = Nu_fd(); xs = np.linspace(0, X, n); dx = xs[1] - xs[0]; Hloc = 2.0 if kind == "step" else 1.0; Uloc = 1.0 / Hloc
    # upstream stretch x_c .. 0 is added by the caller; here start at x = 0 with the given bulk temperature
    a = nu * (1.0 / PE) / (Uloc * Hloc ** 2); Tb = Tb_c * np.exp(-a * xs); q = (nu / Hloc) * Tb      # q in units alpha/H ... times (1/H) already for depth Hloc: q = alpha Nu Tb / Hloc  [alpha units]
    Q = np.concatenate([[0], np.cumsum(0.5 * (q[1:] + q[:-1]) * dx)]); return xs, q, Q
if __name__ == "__main__":
    print(Nu_fd(), PE)
