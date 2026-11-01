"""Restart-safe job queue: `python -m ice_cliff_survival.reactivation.queue OUTDIR JOBFILE NWORKERS` where JOBFILE lines are `name layout dx dtm_h days removal [k=v ...]`. Skips finished results; admits a job only if free memory (MemAvailable) exceeds 1.5 GB; logs wall times."""
import os, subprocess, sys, time, json
def memavail_mb():
    for l in open("/proc/meminfo"):
        if l.startswith("MemAvailable"): return int(l.split()[1]) / 1024
def main(out, jobfile, nw):
    jobs = [l.split() for l in open(jobfile) if l.strip() and not l.startswith("#")]; running = []; t0 = time.time(); log = open(os.path.join(out, "queue.log"), "a"); os.makedirs(out, exist_ok=True)
    while jobs or running:
        running = [(p, n, t) for p, n, t in running if p.poll() is None or log.write(f"{time.time()-t0:8.0f}s done {n} in {time.time()-t:.0f}s rc={p.returncode}\n") is None or log.flush()]
        while jobs and len(running) < nw and memavail_mb() > 1500:
            j = jobs.pop(0)
            if os.path.exists(os.path.join(out, j[0] + ".json")): continue
            p = subprocess.Popen([sys.executable, "-m", "ice_cliff_survival.reactivation.run_case", out, j[0], j[1], j[2], j[3], j[4], j[5]] + j[6:], stdout=open(os.path.join(out, j[0] + ".log"), "w"), stderr=subprocess.STDOUT); running.append((p, j[0], time.time())); log.write(f"{time.time()-t0:8.0f}s start {j[0]}\n"); log.flush()
        time.sleep(2)
    log.write(f"{time.time()-t0:8.0f}s ALL DONE (campaign wall time)\n"); log.close()
if __name__ == "__main__": main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
