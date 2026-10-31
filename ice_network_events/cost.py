"""Cost of each stage (wall s, single process) for one test configuration, plus network solver times from the records. -> cost.json"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import design, sim, methods as M, series as S
c, dev, test = design.selection(); cfg = test["E04"]; p = dict(cfg["geom"], deficit=-0.10); ser = S.make_series(cfg["forcing"], 36.0); out = {}
t = time.perf_counter(); th = M.thresholds(p); out["tabulate_thermal_2_tables_s"] = time.perf_counter() - t
t = time.perf_counter(); ft = sim.fixed_tabs(p); out["tabulate_fixed_2_tables_s"] = time.perf_counter() - t
t = time.perf_counter(); crit = M.criterion(ser, th, p); out["criterion_A_s"] = time.perf_counter() - t
t = time.perf_counter(); r, _ = sim.se_run(p, ser, ft); out["SE_run_s"] = time.perf_counter() - t
se = dict(t=r["t"], wet=r["Q2"] > 0.005 * S.means(ser)["Q_bar"])
for m in ("TR", "AVG", "PW6", "PW3", "EA"):
    t = time.perf_counter(); sim.rtm_run(p, ser, m, th, ft, se); out[f"RTM_{m}_s"] = time.perf_counter() - t
ev = json.load(open(os.path.join(HERE, "eval_test.json")))
out["network_r17_solver_s_median"] = sorted(v["solver_s"] for r_ in ev.values() for k, v in r_["net"].items())[len(ev) * 3 - 1]
json.dump(out, open(os.path.join(HERE, "cost.json"), "w"), indent=1); print(json.dumps(out, indent=1))
