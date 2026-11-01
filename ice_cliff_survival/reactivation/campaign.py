"""Evaluation campaign launcher. `python -m ice_cliff_survival.reactivation.campaign SPEC_DIR OUT_DIR [workers]`. Refuses to run unless (i) production_spec.json matches its recorded sha256, (ii) every solver file hash in the spec matches the working tree, (iii) the reactivation directory has no uncommitted changes. Medium runs (longest first) then coarse repeats for numerical-error checks."""
import hashlib, json, os, subprocess, sys
from . import queue as Q
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def main(spec_dir, out, nw):
    p = os.path.join(spec_dir, "production_spec.json")
    if sha(p) != open(os.path.join(spec_dir, "production_spec.sha256")).read().strip(): sys.exit("REFUSED: spec hash mismatch")
    spec = json.load(open(p))
    for f, h in spec["solver"]["file_sha256"].items():
        if sha(os.path.join(HERE, f)) != h: sys.exit(f"REFUSED: {f} changed since the freeze")
    dirty = subprocess.run(["git", "-C", ROOT, "status", "--porcelain", "ice_cliff_survival/reactivation"], capture_output=True, text=True).stdout.strip()
    if any(l for l in dirty.splitlines() if "results/" not in l and "figures" not in l): sys.exit("REFUSED: uncommitted changes:\n" + dirty)
    os.makedirs(out, exist_ok=True); jobs = []
    for r in spec["evaluation_runs"]: jobs.append(f"{r['name']}_medium {r['layout']} {r['dx']} {r['dtm_h']} {r['days']} {'none' if r['removal_day'] is None else r['removal_day']}")
    for r in spec["evaluation_runs"]: jobs.append(f"{r['name']}_coarse {r['layout']} 0.05 {r['dtm_h']} {r['days']} {'none' if r['removal_day'] is None else r['removal_day']}")
    jf = os.path.join(out, "jobs.txt"); open(jf, "w").write("\n".join(jobs) + "\n"); json.dump(dict(spec_sha256=sha(p), git_head=subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()), open(os.path.join(out, "launch_record.json"), "w"))
    Q.main(out, jf, nw)
if __name__ == "__main__": main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 4)
