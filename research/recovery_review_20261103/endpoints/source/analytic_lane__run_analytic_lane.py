"""Run PROTOCOL_FROZEN.json: ODE traces of the ideal two-route model scored by the packet scorer. python3 run_analytic_lane.py OUT.json"""
import hashlib, importlib.util, json, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
spec_ = importlib.util.spec_from_file_location("score_trace", os.path.join(ROOT, "research", "glacier_campaign_20261103", "scoring", "score_trace.py")); st = importlib.util.module_from_spec(spec_); spec_.loader.exec_module(st)
RHOI, L, RHOW, CP, H, T, Q0, GAP = 917.0, 334000.0, 1000.0, 4217.0, 200.0, 0.5, 0.07, 0.25
K = (2 / 3) ** 1.5 * math.sqrt(9.81) * 0.5; HC = (Q0 / K) ** (2 / 3)
def melt(q, A):
    if q <= 0: return 0.0
    f = RHOW * CP * q; return f * T * (-math.expm1(-H * A / f)) / (RHOI * L * A)
def split(Q, g):
    lo, hi = 0.0, g + (Q / K) ** (2 / 3)
    for _ in range(100):
        e = 0.5 * (lo + hi)
        if K * max(e - g, 0) ** 1.5 + K * e ** 1.5 < Q: lo = e
        else: hi = e
    e = 0.5 * (lo + hi); return K * max(e - g, 0) ** 1.5, K * e ** 1.5
def trace(At, mult, dur_h, dt=900.0, evolve=True):
    ps, pe = 86400.0, 86400.0 + dur_h * 3600.0; end = pe + 172800.0; g = GAP; t = 0.0; rows = []
    Qf = lambda tt: Q0 * mult if ps <= tt < pe else Q0
    def rec(g, t):
        Q = Qf(t); qt, qb = split(Q, max(g, 0.0)); rows.append(dict(time_s=t, toe_fraction=min(max(qt / Q, 0.0), 1.0), delivered_heat_W=RHOI * L * melt(qt, At) * At))
    rate = lambda gg, tt: -(melt(split(Qf(tt), max(gg, 0.0))[0], At) - melt(split(Qf(tt), max(gg, 0.0))[1], 1.0)) if evolve else 0.0
    rec(g, t)
    while t < end - 1e-9:
        k1 = rate(g, t); k2 = rate(g + 0.5 * dt * k1, t + 0.5 * dt); g += dt * k2; t += dt; rec(g, t)
    return rows, g, (ps, pe, end)
def main(out):
    spec_json = json.load(open(os.path.join(HERE, "PROTOCOL_FROZEN.json"))); sha = hashlib.sha256(open(os.path.join(HERE, "PROTOCOL_FROZEN.json"), "rb").read()).hexdigest(); res = dict(protocol_sha256=sha, cases=[], violations=[])
    cap = lambda mult: (H * T / (RHOI * L) - melt(Q0 * mult / 2.0, 1.0)) * 21 * 86400
    for At in (1.0, 0.5, 0.25):
        for mult in (1, 2, 4, 8):
            for dur in (6, 24, 72):
                for ctl in ("evolving", "fixed_geometry"):
                    rows, gend, (ps, pe, end) = trace(At, mult, dur, evolve=(ctl == "evolving")); rows2, gend2, _ = trace(At, mult, dur, dt=450.0, evolve=(ctl == "evolving"))
                    sc = st.score(rows, dict(pre_start_s=0.0, pulse_start_s=ps, pulse_end_s=pe, post_end_s=end, contact_threshold=1e-3, maximum_sampling_gap_s=900.0))
                    peak_rewet = max(r["toe_fraction"] for r in rows if ps <= r["time_s"] <= pe); case = dict(id=f"A{At:g}_Q{mult}_D{dur}h_{ctl}", toe_area_m2=At, pulse_multiple=mult, pulse_h=dur, control=ctl, classification=sc["classification"], peak_toe_fraction_in_pulse=peak_rewet,
                        sustained_postpulse=sc["sustained_over_entire_postpulse_window"], postpulse_wet_s=sc["postpulse"]["wet_seconds"], postpulse_dry_s=sc["postpulse"]["dry_seconds"], postpulse_signed_heat_J=sc["postpulse"]["signed_delivered_heat_J"], gap_end_m=gend, gap_change_m=gend - GAP, gap_end_dt450_m=gend2, dt_difference_m=gend - gend2, analytic_21d_repair_bound_m=cap(mult), prefix_sha256=sc.get("preloss_prefix_sha256") or sc.get("fingerprint") or None)
                    case["_scorer_keys"] = sorted(sc.keys()) if not res["cases"] else None; res["cases"].append(case)
                    if ctl == "evolving":
                        if sc["classification"] == "sustained_postpulse_recovery": res["violations"].append(("sustained", case["id"]))
                        if GAP - gend > cap(mult) + 1e-6: res["violations"].append(("exceeds_bound", case["id"]))
    ev = [c for c in res["cases"] if c["control"] == "evolving"]; fx = [c for c in res["cases"] if c["control"] == "fixed_geometry"]
    res["summary"] = dict(n_evolving=len(ev), classes_evolving={k: sum(c["classification"] == k for c in ev) for k in sorted({c["classification"] for c in ev})}, classes_fixed={k: sum(c["classification"] == k for c in fx) for k in sorted({c["classification"] for c in fx})},
                          max_gap_decrease_m=max(GAP - c["gap_end_m"] for c in ev), min_gap_change_m=min(c["gap_change_m"] for c in ev), cases_with_rewetting_in_pulse=sum(c["peak_toe_fraction_in_pulse"] >= 1e-3 for c in ev),
                          evolving_vs_fixed_max_abs_gap_difference_m=max(abs(a["gap_end_m"] - b["gap_end_m"]) for a, b in zip(ev, fx)), max_abs_dt_difference_m=max(abs(c["dt_difference_m"]) for c in ev), violations=res["violations"])
    json.dump(res, open(out, "w"), indent=1, allow_nan=False); print(json.dumps(res["summary"], indent=1)); print(res["cases"][0]["_scorer_keys"])
if __name__ == "__main__": main(sys.argv[1])
