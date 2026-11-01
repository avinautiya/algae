"""python -m ice_cliff_survival.reactivation.campaign2 SPEC_DIR OUT_DIR [workers]: runs the v2 job list from production_spec_v2.json (refuses on hash mismatch / dirty tree); medium runs first."""
import hashlib, json, os, subprocess, sys
from . import queue as Q
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def main(spec_dir, out, nw):
    p = os.path.join(spec_dir, "production_spec_v2.json")
    if sha(p) != open(os.path.join(spec_dir, "production_spec_v2.sha256")).read().strip(): sys.exit("REFUSED: spec hash mismatch")
    spec = json.load(open(p))
    for f, h in spec["solver"]["file_sha256"].items():
        if f != "campaign2.py" and sha(os.path.join(HERE, f)) != h: sys.exit(f"REFUSED: {f} changed since the freeze")
    dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "ice_cliff_survival/reactivation"], capture_output=True, text=True).stdout
    if any(l for l in dirty.splitlines() if "results/" not in l and "figures" not in l): sys.exit("REFUSED: uncommitted changes\n" + dirty)
    os.makedirs(out, exist_ok=True); runs = sorted(spec["runs"], key=lambda r: (r["dx"], -r["days"]))
    jobs = [f"{r['name']} {r['layout']} {r['dx']} {r['dtm_h']} {r['days']} {'none' if r['removal_day'] is None else r['removal_day']} " + " ".join(f"{k}={v}" for k, v in r["extra"].items()) for r in runs]
    jobs.sort(key=lambda l: float(l.split()[2]))      # medium (0.025) first
    jf = os.path.join(out, "jobs.txt"); open(jf, "w").write("\n".join(jobs) + "\n"); json.dump(dict(spec_sha256=sha(p), git_head=subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()), open(os.path.join(out, "launch_record.json"), "w")); Q.main(out, jf, nw)
if __name__ == "__main__": main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4)
