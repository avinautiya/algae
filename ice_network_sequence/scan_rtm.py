"""RTM_seq scan over deterministic random plausible configurations (seed 5) used ONLY to stratify the test set by predicted effect size; no simulator is run. -> scan_rtm.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sdesign as design, predict_seq
rng = np.random.default_rng(5); out = []
GEOMS = [dict(), dict(S2=0.01), dict(w2=0.4, L2=100.0, S2=0.005), dict(L2=120.0), dict(w2=0.8), dict(S1=0.02, L1=60.0, S2=0.005, L2=100.0)]
for i in range(120):
    g = dict(GEOMS[rng.integers(len(GEOMS))]); g["offset"] = float(rng.choice([0.02, 0.05, 0.08, 0.12, 0.16, 0.20])) if "S1" not in g else 0.05; g["area"] = float(np.exp(rng.uniform(np.log(300), np.log(8000))))
    QA = float(rng.uniform(0.02, 0.06)); TA = float(rng.uniform(1.0, 3.0)); QB = float(rng.uniform(0.06, 0.12)); TB = float(rng.uniform(0.1, 0.6))
    if rng.random() < 0.3: QA, QB = QB, QA          # hot-high vs cold-low variants
    if rng.random() < 0.25: QB = QA                  # temperature-only contrast
    if rng.random() < 0.2: TB = TA                   # discharge-only contrast
    tau = float(rng.choice([3.0, 6.0, 12.0])); c = design.cfg(g, ((QA, TA), (QB, TB)), tau_h=tau)
    try:
        pr = predict_seq.predict(f"X{i:03d}", c); th, fx = pr["thermal"], pr["fixed"]
        out.append(dict(i=i, config=c, share_eff=th["AB"]["share2_end"] - th["BA"]["share2_end"], share_eff_fixed=fx["AB"]["share2_end"] - fx["BA"]["share2_end"], eta=(th["AB"]["heat1"] + th["AB"]["heat2"]) / (th["BA"]["heat1"] + th["BA"]["heat2"]) - 1,
                        dom=[th[k]["dominant_end"] for k in ("AB", "BA", "AVG")], share_avg=th["AVG"]["share2_end"], tie_dist=min(abs(th[k]["share2_end"] - 0.5) for k in ("AB", "BA", "AVG")), sill_eff=[th["AB"]["sill1_mm"] - th["BA"]["sill1_mm"], th["AB"]["sill2_mm"] - th["BA"]["sill2_mm"]]))
    except Exception as e: out.append(dict(i=i, config=c, error=str(e)))
    if i % 10 == 9: json.dump(out, open(os.path.join(HERE, "scan_rtm.json"), "w"), default=float); print(i + 1, flush=True)
json.dump(out, open(os.path.join(HERE, "scan_rtm.json"), "w"), default=float)
