"""Restricted scalar falsifier model:  R = F(T, s),  ds/dt = S0 + D F(T, s) - C(s, Q),  D > 0 constant, Q fixed.
Hypotheses (H): F_T > 0, F_s < 0, C_s >= 0, C(0, Q) = 0 (so s >= 0 is invariant without clipping), F, C continuously differentiable.
Results (proved in MODEL_FOUNDATIONS.md §6, checked numerically in tests): (a) g(s) = S0 + D F - C is strictly decreasing in s, so an interior equilibrium is unique on a connected admissible interval and linearly stable; existence must be checked (bracket).
(b) ds*/dT = D F_T / (C_s - D F_s) > 0 ; dR*/dT = F_T C_s / (C_s - D F_s) >= 0.  (c) for T_H(t) >= T_L(t) pointwise, same s(0), Q(t), S0(t): s_H(t) >= s_L(t) and D int (R_H - R_L) dt = (s_H - s_L) + int (C(s_H) - C(s_L)) dt >= 0.
The integrator below never clips s; it raises if s leaves the admissible interval (violating hypotheses) so a clipped trajectory can never masquerade as a theorem instance."""
from dataclasses import dataclass
from typing import Callable
import math
@dataclass(frozen=True)
class ScalarModel:
    F: Callable; C: Callable; D: float; S0: Callable      # F(T,s), C(s,Q), S0(t)
    dF_dT: Callable = None; dF_ds: Callable = None; dC_ds: Callable = None
    def rhs(self, t, s, T, Q): return self.S0(t) + self.D * self.F(T, s) - self.C(s, Q)
def _d(f, x, h=1e-6): return (f(x + h) - f(x - h)) / (2 * h)
def equilibrium(m, T, Q, s_lo=0.0, s_hi=1e6, t=0.0, tol=1e-14):
    """unique root of g(s)=S0+D F-C on [s_lo,s_hi] if the bracket exists, else None (existence checked separately)."""
    g = lambda s: m.rhs(t, s, T, Q)
    if g(s_lo) < 0 or g(s_hi) > 0: return None
    lo, hi = s_lo, s_hi
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if g(mid) > 0: lo = mid
        else: hi = mid
        if hi - lo < tol * max(1.0, abs(hi)): break
    return 0.5 * (lo + hi)
def sensitivities(m, T, Q, s_star, t=0.0):
    F_T = m.dF_dT(T, s_star) if m.dF_dT else _d(lambda x: m.F(x, s_star), T); F_s = m.dF_ds(T, s_star) if m.dF_ds else _d(lambda x: m.F(T, x), s_star); C_s = m.dC_ds(s_star, Q) if m.dC_ds else _d(lambda x: m.C(x, Q), s_star)
    den = C_s - m.D * F_s; return dict(F_T=F_T, F_s=F_s, C_s=C_s, ds_dT=m.D * F_T / den, dR_dT=F_T * C_s / den, g_s=-den)
def integrate(m, T_of_t, Q_of_t, s0, t_end, n=4000, s_max=1e9):
    """RK4 on the augmented state (s, I_R = int R dt, I_C = int C dt); returns final state and arrays. No clipping: raises if s < 0 or non-finite."""
    h = t_end / n; s = s0; IR = 0.0; IC = 0.0; ts = [0.0]; ss = [s0]
    def f(t, y):
        s_, _, _ = y; T = T_of_t(t); Q = Q_of_t(t); R = m.F(T, s_); Cv = m.C(s_, Q); return (m.S0(t) + m.D * R - Cv, R, Cv)
    y = (s, IR, IC)
    for i in range(n):
        t = i * h; k1 = f(t, y); k2 = f(t + h / 2, tuple(a + h / 2 * b for a, b in zip(y, k1))); k3 = f(t + h / 2, tuple(a + h / 2 * b for a, b in zip(y, k2))); k4 = f(t + h, tuple(a + h * b for a, b in zip(y, k3)))
        y = tuple(a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4))
        if not math.isfinite(y[0]) or y[0] < -1e-9 or y[0] > s_max: raise ValueError(f"trajectory left the admissible interval at t={t+h}: s={y[0]} (hypotheses violated; no clipping applied)")
        ts.append(t + h); ss.append(y[0])
    return dict(s=y[0], int_R=y[1], int_C=y[2], t=ts, s_path=ss)
def accumulation_hitting_bound(Vs0, Vcrit, Smin, Cmax):
    """upper hitting-time bound t* <= (Vcrit - Vs0)/(Smin - Cmax) for Smin > Cmax (supply >= Smin, admissible export <= Cmax, persistently). Not a minimum time; not proof of burial."""
    if Smin <= Cmax: return math.inf
    return max(0.0, (Vcrit - Vs0) / (Smin - Cmax))
