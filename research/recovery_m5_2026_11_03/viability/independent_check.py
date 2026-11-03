"""Independent re-derivation of the ideal two-route viability bound (written without importing screen.py).
Checks (i) lemma assumptions numerically, (ii) the bound, (iii) full ODE integration of the declared ideal model over 21 d for the 12 declared cases
AND for an adversarial sweep of toe area (parameters of the declared model only; nothing is inflated), (iv) the dimensionless dependence of the bound on NTU
(analysis of which assumption controls the result; NOT a proposal to change parameters). Writes viability_check.json.
python3 independent_check.py OUT.json"""
import json, math, sys
RHOI, L, RHOW, CP = 917.0, 334000.0, 1000.0, 4217.0
H, T, Q0, GAP, HORIZON = 200.0, 0.5, 0.07, 0.25, 21.0 * 86400.0
K = (2.0 / 3.0) ** 1.5 * math.sqrt(9.81) * 0.5
HC = (Q0 / K) ** (2.0 / 3.0)
FILM = H * T / (RHOI * L)                                   # m/s, infinite-water limit of the floor melt rate
def melt(q, A, h=H, T_=T):                                  # finite-water Stefan floor lowering rate (m/s) of a route of area A carrying q
    if q <= 0.0: return 0.0
    f = RHOW * CP * q; return f * T_ * (-math.expm1(-h * A / f)) / (RHOI * L * A)
def split(Q, g):                                            # equal-K critical controls; bypass crest at 0, toe crest at g
    lo, hi = 0.0, g + (Q / K) ** (2.0 / 3.0)
    for _ in range(200):
        eta = 0.5 * (lo + hi)
        if K * max(eta - g, 0.0) ** 1.5 + K * eta ** 1.5 < Q: lo = eta
        else: hi = eta
    eta = 0.5 * (lo + hi); return K * max(eta - g, 0.0) ** 1.5, K * eta ** 1.5
def integrate(Q, At, Ab=1.0, g0=GAP, horizon=HORIZON, dt=900.0):
    """RK2 on the gap g=z_toe-z_bypass (m). dg/dt = -(E_toe-E_byp). Returns (g_end, max_repair_so_far, time to reach HC or None)."""
    g, t, best, hit = g0, 0.0, 0.0, None
    def rate(gg):
        qt, qb = split(Q, max(gg, 0.0)); return -(melt(qt, At) - melt(qb, Ab))
    while t < horizon:
        k1 = rate(g); k2 = rate(g + 0.5 * dt * k1); g += dt * k2; t += dt
        best = max(best, g0 - g)
        if hit is None and g <= HC: hit = t
        if g <= 0.0: break
    return g, best, hit
def main(out):
    res = {"assumptions_checked": {}, "cases": [], "adversarial_toe_area_sweep": [], "ntu_dependence": []}
    # (i) lemma checks
    qs = [10 ** (-6 + 0.05 * i) for i in range(160)]
    res["assumptions_checked"]["melt_rate_monotone_increasing_in_q"] = all(melt(q2, 1.0) >= melt(q1, 1.0) for q1, q2 in zip(qs, qs[1:]))
    res["assumptions_checked"]["melt_rate_never_exceeds_film_limit"] = all(melt(q, a) <= FILM * (1 + 1e-12) for q in qs for a in (0.01, 0.25, 1.0, 10.0))
    worst = min(split(Q, g)[1] / Q for Q in (0.07, 0.28, 0.56) for g in [0.0, 1e-3, 0.05, 0.19, 0.25, 1.0, 5.0])
    res["assumptions_checked"]["bypass_share_ge_half_for_all_gaps_ge_0_equal_K"] = worst >= 0.5 - 1e-12
    res["assumptions_checked"]["min_bypass_share_found"] = worst
    res["assumptions_checked"]["film_limit_m_per_day"] = FILM * 86400.0
    res["assumptions_checked"]["critical_head_m"] = HC; res["assumptions_checked"]["K"] = K
    need = GAP - HC
    for At in (1.0, 0.5, 0.25):
        for mult in (1, 2, 4, 8):
            Q = Q0 * mult; cap = FILM - melt(Q / 2.0, 1.0); g_end, best, hit = integrate(Q, At)
            res["cases"].append(dict(id=f"A{At:g}_Q{mult}", bound_m_per_day=cap * 86400, bound_21d_m=cap * HORIZON, required_m=need, bound_over_required=cap * HORIZON / need,
                                      ode_repair_21d_m=best, ode_over_bound=best / (cap * HORIZON) if cap > 0 else None, ode_reaches_critical_gap=hit is not None, ode_final_gap_m=g_end))
    # (iii) adversarial sweep over toe area for the declared Q multiples (no other parameter changed)
    for mult in (1, 2, 4, 8):
        Q = Q0 * mult; cap = FILM - melt(Q / 2.0, 1.0)
        for At in (1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0):
            g_end, best, hit = integrate(Q, At, horizon=HORIZON)
            res["adversarial_toe_area_sweep"].append(dict(Q_mult=mult, toe_area_m2=At, ode_repair_21d_m=best, bound_21d_m=cap * HORIZON, within_bound=best <= cap * HORIZON * (1 + 1e-9), reaches_critical=hit is not None))
    # (iv) dimensionless bound: cap/FILM = 1 - (1-exp(-N))/N with N = h A_b /(rho c Q/2)
    for mult in (1, 2, 4, 8):
        N = H * 1.0 / (RHOW * CP * Q0 * mult / 2.0); ratio = 1.0 - (-math.expm1(-N)) / N
        res["ntu_dependence"].append(dict(Q_mult=mult, NTU_bypass_at_half_flow=N, cap_over_film=ratio, small_N_approx_N_over_2=N / 2.0))
    # NTU that WOULD be needed for a 21 d feasibility (analysis of the assumption, not a proposal): solve cap_over_film = need/(FILM*horizon)
    target = need / (FILM * HORIZON); lo, hi = 1e-9, 50.0
    for _ in range(200):
        m = 0.5 * (lo + hi); r = 1.0 - (-math.expm1(-m)) / m
        if r < target: lo = m
        else: hi = m
    res["required_cap_over_film_for_21d"] = target; res["required_NTU_for_21d"] = 0.5 * (lo + hi)
    res["required_NTU_over_declared_NTU_Q1"] = res["required_NTU_for_21d"] / res["ntu_dependence"][0]["NTU_bypass_at_half_flow"]
    res["all_ode_within_bound"] = all(c["ode_over_bound"] is None or c["ode_over_bound"] <= 1 + 1e-9 for c in res["cases"]) and all(r["within_bound"] for r in res["adversarial_toe_area_sweep"])
    res["any_ode_reaches_critical_gap"] = any(c["ode_reaches_critical_gap"] for c in res["cases"]) or any(r["reaches_critical"] for r in res["adversarial_toe_area_sweep"])
    json.dump(res, open(out, "w"), indent=1, allow_nan=False)
    print(json.dumps({k: res[k] for k in ("all_ode_within_bound", "any_ode_reaches_critical_gap", "required_NTU_for_21d", "required_NTU_over_declared_NTU_Q1")}), res["assumptions_checked"])
    for c in res["cases"][:4]: print(c)
    print(res["ntu_dependence"])
if __name__ == "__main__": main(sys.argv[1])
