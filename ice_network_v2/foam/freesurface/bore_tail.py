"""Finish the bore study after the driver was stopped: analyse n=8 when interFoam has finished; run the time-step check at n=4 (Co 0.15) which is cheap; records n=12 and n=8/Co0.15 as NOT RUN (cost ~1-2 h each on one core)."""
import json, os, sys, shutil, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import bore as B, common as C, jump as J
wd = sys.argv[1]; out = json.load(open(os.path.join(HERE, "bore_results.json")))
case = f"{wd}/bore_n8_co0.3"; g = B.make.__wrapped__(case, 8) if hasattr(B.make, "__wrapped__") else dict(dy=B.H0 / 8, nx=int(round(B.L / (B.H0 / 8))), ny=int(round(B.HT / (B.H0 / 8))))
while not os.path.exists(f"{case}/1.5/alpha.water"): time.sleep(10)
time.sleep(5); out["n8_co0.3"] = B.analyse(case, g, [0.5, 0.75, 1.0, 1.25, 1.5], 1.5); out["n8_co0.3"]["dy"] = g["dy"]
case4 = f"{wd}/bore_n4_co0.15"; g4 = B.make(case4, 4, maxco=0.15); C.run(["blockMesh"], case4, "bm.log"); shutil.copy(f"{case4}/0/alpha.water.orig", f"{case4}/0/alpha.water"); J.write_internal(case4, g4); C.run(["checkMesh"], case4, "cm.log"); C.run(["interFoam"], case4, "run.log")
out["n4_co0.15"] = B.analyse(case4, g4, [0.5, 0.75, 1.0, 1.25, 1.5], 1.5); out["n4_co0.15"]["dy"] = g4["dy"]
out["n12_co0.3"] = dict(not_run="estimated 1-2 h single core (n=8 took ~35 min of CPU for 1.4 s); not run in this milestone"); out["n8_co0.15"] = dict(not_run="same reason")
json.dump(out, open(os.path.join(HERE, "bore_results.json"), "w"), indent=1); print("done")
