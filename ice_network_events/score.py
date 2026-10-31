"""Apply the frozen rules (FREEZE.md) to eval_test.json -> scores_test.json, SCORES_TEST.md.  python score.py [dev|test]"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); which = sys.argv[1] if len(sys.argv) > 1 else "test"
E = json.load(open(os.path.join(HERE, f"eval_{which}.json"))); M = ("AVG", "PW6", "PW3", "EA", "PWn"); L = []; out = {}
sd = lambda n, m: abs(n[m]["share_deliv"] - n["TR"]["share_deliv"]) * 100       # pp
wh = lambda n, m: abs(n[m]["wet_h"] - n["TR"]["wet_h"])
rows = {}
for cid, r in E.items():
    n = r.get("net")
    if not n or any(m not in n for m in ("TR", "AVG")): continue
    pos = bool(n["AVG"]["cls"] != n["TR"]["cls"] or wh(n, "AVG") > 2 or sd(n, "AVG") > 5); unsafe = bool(r["crit"]["unsafe"])
    rows[cid] = dict(category=r["category"], near=r["near"], R=r["crit"]["R"], unsafe=unsafe, positive=pos, cls={m: n[m]["cls"] for m in n}, wet_err_h={m: wh(n, m) for m in n if m != "TR"}, share_err_pp={m: sd(n, m) for m in n if m != "TR"},
                     final_share_err_pp={m: abs(n[m]["share_end"] - n["TR"]["share_end"]) * 100 for m in n if m != "TR"}, rtm_share_err_pp=abs(r["rtm"]["TR"]["share_deliv"] - n["TR"]["share_deliv"]) * 100,
                     rtm_wet_err_h=abs(r["rtm"]["TR"]["wet_h"] - n["TR"]["wet_h"]), rel_vol_err={m: abs(n[m]["cum2"] - n["TR"]["cum2"]) / max(n["TR"]["cum2"], 1e-9) for m in n if m != "TR"} if n["TR"]["cum2"] > 0.01 * (n["TR"]["cum1"] + n["TR"]["cum2"]) else None, rtm_cls_ok=r["rtm"]["TR"]["cls"] == n["TR"]["cls"], solver_s={m: n[m]["solver_s"] for m in n})
tp = sum(v["unsafe"] and v["positive"] for v in rows.values()); fp = sum(v["unsafe"] and not v["positive"] for v in rows.values()); fn = sum((not v["unsafe"]) and v["positive"] for v in rows.values()); tn = sum((not v["unsafe"]) and (not v["positive"]) for v in rows.values())
P = tp + fn; N = fp + tn; rec = tp / P if P else None; prec = tp / (tp + fp) if tp + fp else None
A_ok = bool(P >= 3 and N >= 3 and rec is not None and prec is not None and rec >= 0.8 and prec >= 0.8)
pos = {k: v for k, v in rows.items() if v["positive"]}
med = lambda m: float(np.median([v["share_err_pp"].get(m, np.nan) for v in pos.values()])) if pos else None
mm = {m: med(m) for m in M}; worse = sum(v["share_err_pp"]["EA"] > v["share_err_pp"]["AVG"] for v in pos.values())
ea_ok = bool(pos and mm["EA"] <= 0.5 * mm["AVG"] and mm["EA"] <= 0.8 * np.nanmin([mm["PW6"], mm["PW3"]]) and worse <= 1)
rtm_med = float(np.median([v["rtm_share_err_pp"] for v in pos.values()])) if pos else None
relm = {m: float(np.median([v["rel_vol_err"][m] for v in rows.values() if v["rel_vol_err"] and v["positive"] and m in v["rel_vol_err"]] or [np.nan])) for m in M}
out = dict(median_rel_vol_err_positive=relm, confusion=dict(tp=tp, fp=fp, fn=fn, tn=tn), recall=rec, precision=prec, criterion_A_established=A_ok, median_share_err_pp=mm, EA_worse_than_AVG_count=worse, n_positive=len(pos), EA_useful_improvement=ea_ok, RTM_TR_median_share_err_pp=rtm_med, RTM_TR_success=bool(rtm_med is not None and rtm_med <= 2), rows=rows)
# robustness
rob = {}
for cid, r in E.items():
    for key in ("ref", "net_nu0.5", "net_nu2"):
        d = r.get(key)
        if d and "TR" in d and "AVG" in d:
            ref_net = r["net"]; ok_cls = all(d[m]["cls"] == ref_net[m]["cls"] for m in d if m in ref_net)
            order = (abs(d["EA"]["share_deliv"] - d["TR"]["share_deliv"]) < abs(d["AVG"]["share_deliv"] - d["TR"]["share_deliv"])) if "EA" in d else None
            order0 = abs(ref_net["EA"]["share_deliv"] - ref_net["TR"]["share_deliv"]) < abs(ref_net["AVG"]["share_deliv"] - ref_net["TR"]["share_deliv"])
            rob.setdefault(cid, {})[key] = dict(class_preserved=ok_cls, EA_beats_AVG=order, EA_beats_AVG_r17=order0, share_deliv_pp={m: d[m]["share_deliv"] * 100 for m in d}, share_deliv_r17_pp={m: ref_net[m]["share_deliv"] * 100 for m in d if m in ref_net})
    if r.get("net_legacytol"): rob.setdefault(cid, {})["legacy_tol_TR_share_pp"] = [r["net_legacytol"]["share_deliv"] * 100, r["net"]["TR"]["share_deliv"] * 100]
out["robustness"] = rob; json.dump(out, open(os.path.join(HERE, f"scores_{which}.json"), "w"), indent=1, default=float)
L.append(f"# Scores ({which}) — frozen rules\n"); L.append(f"Criterion A: tp {tp} fp {fp} fn {fn} tn {tn}; recall {rec}, precision {prec}; established: **{A_ok}**\n")
L.append("| cfg | category | R | unsafe | AVG misrepresents | class TR/AVG/PW6/PW3/EA | share err (pp) AVG / PW6 / PW3 / PWn / EA | wet err (h) AVG / PW6 / PW3 / EA | RTM-TR share err |\n|---|---|---|---|---|---|---|---|---|")
for cid, v in rows.items():
    c = v["cls"]; s = v["share_err_pp"]; w = v["wet_err_h"]
    L.append(f"| {cid} | {v['category']}{'*' if v['near'] else ''} | {v['R']:.2f} | {v['unsafe']} | {v['positive']} | {'/'.join(c.get(m,'-')[:3] for m in ('TR','AVG','PW6','PW3','EA'))} | " + " / ".join(f"{s.get(m, float('nan')):.1f}" for m in ('AVG','PW6','PW3','PWn','EA')) + " | " + " / ".join(f"{w.get(m, float('nan')):.1f}" for m in ('AVG','PW6','PW3','EA')) + f" | {v['rtm_share_err_pp']:.2f} |")
L.append(f"\nTruth-positive configs: {len(pos)}. Median |share error| (pp): " + ", ".join(f"{m} {mm[m]:.2f}" for m in M) + f". EA worse than AVG on {worse}. **EA useful improvement (rule 2): {ea_ok}**. RTM-TR median error {rtm_med} pp → success {out['RTM_TR_success']}.")
L.append("\n## Robustness\n```\n" + json.dumps(rob, indent=1, default=float)[:6000] + "\n```")
open(os.path.join(HERE, f"SCORES_{which.upper()}.md"), "w").write("\n".join(L)); print("\n".join(L[:30]))
