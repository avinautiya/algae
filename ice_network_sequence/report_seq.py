"""Generate SCORES_TEST.md and summary statistics from scores_test.json / eval_*.json (all numbers from records)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); S = json.load(open(os.path.join(HERE, "scores_test.json"))); ET = json.load(open(os.path.join(HERE, "eval_test.json"))); ED = json.load(open(os.path.join(HERE, "eval_dev.json")))
L = []; P = L.append
f = lambda o, k, n=4: "-" if (not o or o.get(k) is None) else f"{o[k]:+.{n}f}"
P("## Frozen test: sequence effect AB - BA (predicted by RTM_seq before the runs; observed in the network and the resolved reference)\n")
P("| id | config | pred share_end | net r17 | net r33 | ref 200 | ref 400 | fixed-geometry net | pred heat rel | net heat rel | ref heat rel | pred meaningful routing | observed (net / ref) | dominant AB/BA/AVG (net) | tie | flip |"); P("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for cid, r in S.items():
    n, ref = r["net17"], r["ref400"] or r["ref200"]
    P(f"| {cid} | {r['description'].split(':')[0][:16]} | {r['predicted']['share_end']:+.4f} | {f(n,'share_end')} | {f(r['net33'],'share_end')} | {f(r['ref200'],'share_end')} | {f(r['ref400'],'share_end')} | {f(r['net17_fixed'],'share_end')} | {r['predicted']['heat_rel']:+.3f} | {f(n,'heat_rel',3)} | {f(ref,'heat_rel',3)} | {r['predicted']['meaningful_routing']} | {n['meaningful_routing'] if n else '-'} / {ref['meaningful_routing'] if ref else '-'} | {'/'.join(map(str, n['dominant'])) if n else '-'} | {n['tie'] if n else '-'} | {n['flip'] if n else '-'} |")
P("\n## Cumulative delivered-share effect, persistence, incision, averaging error (network r17)\n"); P("| id | delivered-share effect | cum2 effect (m3) | persistence (end/start-of-continuation share effect) | incision rel | pred avg error (AVG-AB, AVG-BA) | observed avg error | averaging adequate | |AVG - mean(AB,BA)| / |AB - BA| (share_end) |"); P("|---|---|---|---|---|---|---|---|---|")
for cid, r in S.items():
    n = r["net17"]; ab = ET[cid]["thermal_r17"]; mid = 0.5 * (ab["AB"]["share_end"] + ab["BA"]["share_end"]); ratio = abs(ab["AVG"]["share_end"] - mid) / max(abs(n["share_end"]), 1e-4)
    P(f"| {cid} | {n['share_deliv']:+.4f} | {n['cum2']:+.0f} | {'-' if n['persistence'] is None else '%.2f' % n['persistence']} | {n['incision_rel']:.3f} | {[round(x, 3) for x in r['predicted']['avg_error_share_end']]} | {[round(x, 3) for x in n['avg_error_share_end']]} | {n['avg_adequate']} | {ratio:.1f} |")
P("\n## Closure multipliers and tolerance (network r17)\n"); P("| id | base share effect | nu x0.5 | nu x2 | base heat rel | nu x0.5 | nu x2 | destination flip (any) | legacy vs corrected tolerance (share effect) |"); P("|---|---|---|---|---|---|---|---|---|")
for cid in ("X04", "X08", "X09"):
    r = S[cid]; c1, c2 = r["closure"]["nu0.5"], r["closure"]["nu2"]; t = r["tol"]
    P(f"| {cid} | {f(r['net17'],'share_end')} | {f(c1,'share_end')} | {f(c2,'share_end')} | {f(r['net17'],'heat_rel',3)} | {f(c1,'heat_rel',3)} | {f(c2,'heat_rel',3)} | {any(x['flip'] for x in (c1, c2) if x)} | {f(t,'share_end')} vs {f(r['net17'],'share_end')} |" if t is None or True else "")
P("\n## Numerical uncertainty of the effect (3 x rule)\n")
for cid, r in S.items():
    if r["numerical_uncertainty"]: P(f"* {cid}: " + ", ".join(f"{k} {v:.4f}" for k, v in r["numerical_uncertainty"].items()))
# summary statistics
res = [(cid, r["predicted"]["share_end"], r["net17"]["share_end"], (r["ref400"] or r["ref200"] or {}).get("share_end")) for cid, r in S.items()]
res_sign = [(c, p, o, rf) for c, p, o, rf in res if abs(o) > 0.003]
sign_ok = sum(np.sign(p) == np.sign(o) for c, p, o, rf in res_sign)
ratios = [p / o for c, p, o, rf in res_sign]; ratios_ref = [p / rf for c, p, o, rf in res_sign if rf is not None and abs(rf) > 0.003]
P(f"\n## Summary\n* sign of the share effect agrees in {sign_ok} of {len(res_sign)} configurations with |observed effect| > 0.003 (network); against the resolved reference {sum(np.sign(p) == np.sign(rf) for c, p, o, rf in res_sign if rf is not None)} of {sum(rf is not None for c, p, o, rf in res_sign)}")
P(f"* predicted/observed magnitude (network): median {np.median(ratios):.2f}, range {min(ratios):.2f}-{max(ratios):.2f}; (reference): median {np.median(ratios_ref):.2f}, range {min(ratios_ref):.2f}-{max(ratios_ref):.2f}")
hp = [(r["predicted"]["heat_rel"], r["net17"]["heat_rel"]) for r in S.values() if r["net17"]["heat_rel"] is not None]; he = [abs(a - b) for a, b in hp]; P(f"* heat-difference prediction (RTM_seq) vs network: max |error| {max(he):.3f}, median {np.median(he):.4f}")
et = []
for cid, r in S.items():
    pass
P("* destination flips observed (network r17): " + str(sum(r['net17']['flip'] for r in S.values())) + " of 10; resolved reference (informative cases): " + str(sum(bool((r['ref400'] or r['ref200'] or {}).get('flip')) for r in S.values())))
P("* meaningful routing effect predicted for " + str([c for c, r in S.items() if r['predicted']['meaningful_routing']]) + ", observed (network) for " + str([c for c, r in S.items() if r['net17']['meaningful_routing']]) + ", reference " + str([c for c, r in S.items() if (r['ref400'] or r['ref200'] or {}).get('meaningful_routing')]))
P("* meaningful incision effect predicted for " + str([c for c, r in S.items() if r['predicted']['meaningful_incision']]) + ", observed (network) for " + str([c for c, r in S.items() if r['net17']['meaningful_incision']]))
open(os.path.join(HERE, "SCORES_TEST.md"), "w").write("\n".join(L)); print("\n".join(L[-14:]))
