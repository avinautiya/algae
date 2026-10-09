#!/usr/bin/env python
"""
Dependency-aware Phase 1 job registry and runner (replaces the single blocking shell queue).

    python phase1/jobs.py status                 # table of jobs, states, validation
    python phase1/jobs.py run [--max-jobs 2]     # launch every runnable job; returns when none can start
                                                 # and all launched jobs have exited
    python phase1/jobs.py run --only L2_B3LYP_TDA15
    python phase1/jobs.py show L2_CAM_TDA15      # exact command, outdir, dependencies

Durability model (docs/phase1_jobs.md):
  * every job writes into results/v2/<id>/ (persistent disk); excited-state solves checkpoint into
    <outdir>/td_ckpt/ (phase1/tdcheckpoint.py), so a killed job resumes from its last chunk;
  * job state lives in results/v2/jobs_state.json (atomic writes); a job recorded as running whose
    PID is gone is marked 'interrupted' and relaunched on the next `run`;
  * a job is 'complete' only when completion.validate_run accepts its outdir (COMPLETE.json with
    matching checksums; production status, or provisional_geometry where the job allows it); any
    other exit is 'failed' with the reasons recorded - the runner's exit code is non-zero then;
  * concurrency is limited by CPU threads and by MemAvailable (/proc/meminfo) including processes
    the runner did not start (e.g. a legacy run).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common"))
import resources as RS  # noqa: E402
import completion  # noqa: E402

R = os.path.join(HERE, "results")
V2 = os.path.join(R, "v2")
STATE = os.path.join(V2, "jobs_state.json")
L2_GEOM = os.path.join(R, "level2", "opt_B3LYP_final.xyz")
LEGACY_CAM_GUESS = os.path.join(R, "legacy_2026-10-09", "level2_cam", "tda_CAM-B3LYP_x0.npz")


def _l2_td(func, nst, tda=True, geom=L2_GEOM, extra=()):
    a = ["--level", "level2", "--skip-opt", "--start-xyz", geom, "--tddft-functionals", func,
         "--nstates", str(nst), "--td-conv-tol", "1e-5"]
    return a + (["--tda"] if tda else []) + list(extra)


JOBS = {
    # production reference (legacy full TD, 30 roots): audited in place, not recomputed
    "L2_B3LYP_FULL30": dict(kind="legacy_audit", outdir=os.path.join(R, "level2"), deps=[], threads=0, mem_mb=0,
                            allow_provisional=True),
    # functional check (B: B3LYP-TDA vs CAM-TDA) and TDA check (A: B3LYP full vs B3LYP TDA), same geometry
    "L2_CAM_TDA15": dict(args=_l2_td("CAM-B3LYP", 15, extra=["--import-guess", LEGACY_CAM_GUESS]),
                         deps=[], threads=2, mem_mb=7000, allow_provisional=True),
    "L2_B3LYP_TDA15": dict(args=_l2_td("B3LYP", 15), deps=[], threads=2, mem_mb=7000, allow_provisional=True),
    # root-count check for the interpreted interval
    "L2_B3LYP_TDA25": dict(args=_l2_td("B3LYP", 25), deps=["L2_B3LYP_TDA15"], threads=2, mem_mb=7000,
                           allow_provisional=True),
    # geometry: continue the Level 2 optimisation, then repeat the TDA spectrum at the relaxed geometry
    "L2_OPT": dict(args=["--level", "level2", "--start-xyz", L2_GEOM, "--opt-functional", "B3LYP", "--no-td",
                         "--maxsteps", "60"], deps=[], threads=2, mem_mb=7000, allow_provisional=False,
                   accept_no_td=True),
    "L2_B3LYP_TDA15_RELAXED": dict(args=None, deps=["L2_OPT"], threads=2, mem_mb=7000, allow_provisional=False,
                                   geometry_from="L2_OPT"),
    # Level 1 core: optimisation + full TD with both functionals (small molecule)
    "L1_FULL": dict(args=["--level", "level1", "--tddft-functionals", "B3LYP", "CAM-B3LYP", "--nstates", "30",
                          "--td-conv-tol", "1e-5"], deps=[], threads=2, mem_mb=3000, allow_provisional=False),
    # Fe(III)-purpurogallin (exploratory mechanistic investigation; see docs/repair_ledger.md P1-FE-*)
    "FE_CAT": dict(args=["--level", "level3_catecholate", "--tddft-functionals", "B3LYP", "--tda", "--nstates", "40",
                         "--td-conv-tol", "1e-5"], deps=[], threads=2, mem_mb=6000, allow_provisional=False),
    # protonation-state sensitivity (carboxylate anion), compared with L2_B3LYP_TDA15 by compare_states.py
    "L2_COO_OPT_TDA15": dict(args=["--level", "level2_carboxylate", "--start-xyz",
                                   os.path.join(HERE, "results", "v2", "inputs", "level2_carboxylate_start.xyz"),
                                   "--opt-functional", "B3LYP", "--maxsteps", "60", "--tddft-functionals", "B3LYP",
                                   "--tda", "--nstates", "15", "--td-conv-tol", "1e-5"],
                             deps=[], threads=2, mem_mb=7000, allow_provisional=True,
                             exploratory="protonation sensitivity; geometry may stop unconverged at 60 steps"),
    "FE_TROP_XTB": dict(args=["--level", "level3_tropolonate", "--xtb-geometry", "--tddft-functionals", "B3LYP",
                              "--tda", "--nstates", "40", "--td-conv-tol", "1e-5"], deps=[], threads=2, mem_mb=6000,
                        allow_provisional=True, exploratory="xTB geometry only"),
}


# Launch order (docs/chemistry_sensitivity_priorities.md): sensitivities the calibration cannot absorb (Fe complex,
# protonation) before checks it absorbs (TDA, root count, geometry). Definitions are unchanged.
PRIORITY = ["L2_B3LYP_FULL30", "L2_CAM_TDA15", "L1_FULL", "FE_CAT", "L2_COO_OPT_TDA15", "FE_TROP_XTB",
            "L2_OPT", "L2_B3LYP_TDA15_RELAXED", "L2_B3LYP_TDA15", "L2_B3LYP_TDA25"]
JOBS = {k: JOBS[k] for k in PRIORITY + [k for k in JOBS if k not in PRIORITY]}

def outdir(jid):
    return JOBS[jid].get("outdir") or os.path.join(V2, jid)


def command(jid):
    j = JOBS[jid]
    if j.get("cmd"):                      # generic command (used by tests)
        return list(j["cmd"])
    args = j["args"]
    if j.get("geometry_from"):
        src = outdir(j["geometry_from"])
        args = _l2_td("B3LYP", 15, geom=os.path.join(src, "opt_B3LYP_final.xyz"))
    args = list(args)
    conv = converged_geometry(jid, args)
    if conv:
        # the job's own optimisation already converged (record + final geometry): do not repeat it on
        # every relaunch (each repeat cost an SCF + gradient + geomeTRIC check, longer than the time
        # between container reboots). The geometry is accepted through its sidecar (completion.py).
        if "--start-xyz" in args:
            args[args.index("--start-xyz") + 1] = conv
        else:
            args += ["--start-xyz", conv]
        args.append("--skip-opt")
    elif "--skip-opt" not in args and "--xtb-geometry" not in args:
        # optimisation jobs resume from their own latest geometry (history is appended, not truncated)
        func = args[args.index("--opt-functional") + 1] if "--opt-functional" in args else "B3LYP"
        last = os.path.join(outdir(jid), f"opt_{func}_last.xyz")
        if os.path.isfile(last):
            if "--start-xyz" in args:
                args[args.index("--start-xyz") + 1] = last
            else:
                args += ["--start-xyz", last]
    if "--no-td" not in args and "--td-chunk" not in args:
        # checkpoint every 2 Davidson iterations: under CPU contention a 5-iteration chunk of a Level 2
        # solve took > 2 h, so a killed job lost all of it (P1-ENV-3)
        args += ["--td-chunk", str(TD_CHUNK)]
    return [sys.executable, os.path.join(HERE, "run_phase1.py"), *args, "--outdir", outdir(jid),
            "--max-memory", str(int(j["mem_mb"] * 0.7)), "--verbose", "4"]


TD_CHUNK = 2


def converged_geometry(jid, args):
    """Path of the job's converged final geometry if its optimisation record says geomeTRIC converged,
    writing the acceptance sidecar <xyz>.geometry.json from that record (status, thresholds, energy and the
    xyz checksum). None for jobs without an optimisation or with an unconverged one."""
    if "--skip-opt" in args or "--xtb-geometry" in args or "--no-td" in args:
        return None
    func = args[args.index("--opt-functional") + 1] if "--opt-functional" in args else "B3LYP"
    od = outdir(jid)
    rec_p, final = os.path.join(od, f"opt_{func}_opt_record.json"), os.path.join(od, f"opt_{func}_final.xyz")
    if not (os.path.isfile(rec_p) and os.path.isfile(final)):
        return None
    rec = json.load(open(rec_p))
    if not rec.get("converged"):
        return None
    import hashlib
    side = final + ".geometry.json"
    sha = hashlib.sha256(open(final, "rb").read()).hexdigest()
    if not os.path.isfile(side) or json.load(open(side)).get("xyz_sha256") != sha:
        tmp = side + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(dict(status="converged", source=rec_p, optimizer=rec.get("optimizer"),
                           thresholds=rec.get("thresholds"), final_energy_Eh=rec.get("final_energy_Eh"),
                           functional=rec.get("functional"), basis=rec.get("basis"), solvent=rec.get("solvent"),
                           xyz_sha256=sha, note="written by jobs.py from the job's own converged optimisation record"),
                      fh, indent=1)
        os.replace(tmp, side)
    return final


def scratch_dir(jid):
    """Per-job scratch (TMPDIR) for PySCF's density-fitting tensors (~4 GB at Level 2). Kept inside the job
    folder so that files of a killed job are found and removed before it is relaunched, instead of piling
    up in /tmp (orphaned files filled the disk once; P1-ENV-3)."""
    return os.path.join(outdir(jid), "scratch")


def clean_scratch(jid):
    import shutil
    d = scratch_dir(jid)
    if os.path.isdir(d):
        shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d, exist_ok=True)
    return d


# --------------------------------------------------------------------------- state
def load_state():
    try:
        return json.load(open(STATE))
    except (OSError, ValueError):
        return {}


def save_state(st):
    os.makedirs(V2, exist_ok=True)
    tmp = STATE + f".tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        json.dump(st, fh, indent=1)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, STATE)


def pid_alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, TypeError, ValueError):
        return False


def mem_available_mb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 1024
    return 0.0


def mem_total_mb():
    for line in open("/proc/meminfo"):
        if line.startswith("MemTotal:"):
            return int(line.split()[1]) / 1024
    return 0.0


def external_phase1_rss_mb(own_pids):
    tot = 0.0
    for p in os.listdir("/proc"):
        if not p.isdigit() or int(p) in own_pids:
            continue
        try:
            cmd = open(f"/proc/{p}/cmdline", "rb").read().split(b"\0")
            if any(c.endswith(b"run_phase1.py") for c in cmd):
                for line in open(f"/proc/{p}/status"):
                    if line.startswith("VmRSS:"):
                        tot += int(line.split()[1]) / 1024
        except OSError:
            continue
    return tot


def external_phase1_threads(own_pids):
    """Threads requested by run_phase1 processes this runner did not start (e.g. legacy runs)."""
    n = 0
    for p in os.listdir("/proc"):
        if not p.isdigit() or int(p) in own_pids:
            continue
        try:
            cmd = open(f"/proc/{p}/cmdline", "rb").read().split(b"\0")
            env = open(f"/proc/{p}/environ", "rb").read().split(b"\0")
        except OSError:
            continue
        if any(c.endswith(b"run_phase1.py") for c in cmd):
            omp = [e for e in env if e.startswith(b"OMP_NUM_THREADS=")]
            n += int(omp[0].split(b"=")[1]) if omp else os.cpu_count()
    return n


def validated(jid):
    j = JOBS[jid]
    ok, st, why = completion.validate_run(outdir(jid), allow_provisional=j.get("allow_provisional", False))
    if not ok and j.get("accept_no_td") and st in ("production", "provisional_geometry", "missing", "diagnostic_only"):
        # optimisation-only job: complete when a converged final geometry exists
        rec = os.path.join(outdir(jid), "opt_B3LYP_opt_record.json")
        if os.path.isfile(rec) and json.load(open(rec)).get("converged"):
            return True, "converged_geometry", []
        return False, st, why + ["optimisation not converged"]
    return ok, st, why


# --------------------------------------------------------------------------- legacy audit
def legacy_audit(jid):
    """Validate a legacy run in place and write COMPLETE.json/DIAGNOSTIC.json for it."""
    import numpy as np
    import pandas as pd
    od = outdir(jid)
    summ = json.load(open(os.path.join(od, "summary.json")))
    arts, ok, problems = [os.path.join(od, "summary.json")], True, []
    log = open(os.path.join(od, "run.log"), errors="replace").read()
    for func in summ.get("tddft", {}):
        stem = os.path.join(od, f"{summ['level']}_{func}")
        st = pd.read_csv(stem + "_states.csv")
        n = summ["nstates"]
        warn = [l for l in log.splitlines() if "TD-DFT roots not converged" in l and "print(" not in l]
        p = completion.td_problems(st.Energy_eV, st.Oscillator_Strength, [not warn] * len(st), n)
        if warn:
            p.append("legacy log reports unconverged roots")
        problems += [f"{func}: {x}" for x in p]
        ok = ok and not p
        arts += [stem + "_states.csv", stem + "_spectrum.csv", stem + "_spectrum.npz"]
    sc = json.load(open(os.path.join(od, "stop_criterion.json")))
    geom = dict(status=sc.get("geometry_status", "not_converged"), source="stop_criterion.json",
                detail=dict(met=sc.get("met"), not_met=sc.get("not_met")))
    summ.setdefault("td_conv_tol", 1e-6)
    summ["legacy_audit"] = dict(problems=problems, note="legacy run: per-root convergence inferred from the "
                                "absence of the solver's non-convergence warning in run.log")
    st = completion.write_marker(od, summ, arts, ok, geom)
    return 0 if st in ("production", "provisional_geometry") else 4


# --------------------------------------------------------------------------- runner
MAX_ATTEMPTS = 3
# memory kept free for analysis work alongside the chemistry (an L2 CAM-B3LYP job was OOM-killed when two
# 7 GB jobs left ~1.2 GB for everything else)
MEM_RESERVE_MB = 3000


def runnable(jid, st, tried=()):
    rec = st.get(jid, {})
    if rec.get("state") in ("complete", "running") or jid in tried:
        return False
    if rec.get("state") == "failed" and rec.get("attempts", 0) >= MAX_ATTEMPTS:
        return False                     # needs an explicit --retry-failed
    return all(validated(d)[0] for d in JOBS[jid]["deps"])


def admit(jid, j):
    """Reserve the job's declared peak memory, threads and scratch (8 GB at Level 2) in the shared budget."""
    scratch_gb = 8.0 if j["mem_mb"] >= 5000 else 2.0
    ok, why, _ = RS.try_reserve(f"job:{jid}", os.getpid(), j["mem_mb"], j["threads"], scratch_gb, outdir(jid)
                                if os.path.isdir(outdir(jid)) else HERE)
    return ok, why


def job_alive(jid, rec):
    """The recorded job process is still THIS job's process: same boot, same start time and command line
    (identity recorded at launch), and its command line names the job's output directory. PID existence
    alone is never accepted: after a reboot PIDs are reused by unrelated processes."""
    if jid not in JOBS:
        return False
    return RS.identity_ok(rec.get("ident") or dict(pid=rec.get("pid")), started=rec.get("started"),
                          must_contain=outdir(jid))


def register_running(st):
    """Jobs still running from an earlier runner keep a reservation (declared peak) under their own pid."""
    for jid, rec in st.items():
        if rec.get("state") == "running" and job_alive(jid, rec):
            with RS._locked() as led:
                led[f"job:{jid}"] = dict(pid=int(rec["pid"]), mem_mb=float(JOBS[jid]["mem_mb"]),
                                         threads=int(JOBS[jid]["threads"]), scratch_gb=8.0, started=time.time(),
                                         ident=rec.get("ident") or RS.process_identity(rec["pid"]))


def run(max_jobs=2, only=None, poll=30, retry_failed=False):
    st = load_state()
    own, tried, waiting_logged = {}, set(), {}
    register_running(st)
    if retry_failed:
        for rec in st.values():
            if rec.get("state") == "failed":
                rec["attempts"] = 0
    for jid, rec in st.items():                       # interrupted jobs from a previous runner
        if rec.get("state") == "running" and not job_alive(jid, rec):
            rec["state"] = "interrupted"           # dead, or its pid now belongs to another process
    for jid in JOBS:                                   # already valid outputs count as complete
        if validated(jid)[0]:
            st.setdefault(jid, {})["state"] = "complete"
    save_state(st)
    ids = [only] if only else list(JOBS)
    cpu = os.cpu_count() or 4
    while True:
        for jid in ids:
            if len(own) >= max_jobs or not runnable(jid, st, tried):
                continue
            j = JOBS[jid]
            if j.get("kind") == "legacy_audit":
                tried.add(jid)
                rc = legacy_audit(jid)
                ok, s, why = validated(jid)
                st[jid] = dict(state="complete" if ok else "failed", status=s, reasons=why, rc=rc, t=time.time())
                save_state(st)
                continue
            # shared, cgroup-aware budget (common/resources.py): the reservation is made atomically for
            # this runner's pid and moved to the job's pid once it is started
            ok_b, why_b = admit(jid, j)
            if not ok_b:
                if waiting_logged.get(jid) != why_b:
                    print(f"[jobs] {jid}: waiting for resources ({'; '.join(why_b)})", flush=True)
                    waiting_logged[jid] = why_b
                continue
            od = outdir(jid)
            os.makedirs(od, exist_ok=True)
            env = RS.thread_env(j["threads"])
            env.update(TMPDIR=clean_scratch(jid), PYSCF_TMPDIR=scratch_dir(jid))
            logf = open(os.path.join(od, "job.log"), "a")
            logf.write(f"\n=== {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(command(jid))}\n")
            logf.flush()
            pr = subprocess.Popen(command(jid), stdout=logf, stderr=subprocess.STDOUT, env=env,
                                  cwd=HERE, start_new_session=True)
            own[jid] = pr
            RS.update_pid(f"job:{jid}", pr.pid)
            tried.add(jid)
            attempts = st.get(jid, {}).get("attempts", 0) + 1
            st[jid] = dict(state="running", pid=pr.pid, started=time.time(), cmd=command(jid), attempts=attempts,
                           ident=RS.process_identity(pr.pid), outdir=od)
            save_state(st)
        done = [k for k, q in own.items() if q.poll() is not None]
        for k in done:
            rc = own.pop(k).returncode
            RS.release(f"job:{k}")
            ok, s, why = validated(k)
            st[k] = dict(state="complete" if ok else "failed", status=s, reasons=why, rc=rc, ended=time.time(),
                         attempts=st.get(k, {}).get("attempts", 1))
            save_state(st)
        if not own and not any(runnable(j, st, tried) for j in ids):
            break
        time.sleep(poll)
    failed = [k for k in ids if st.get(k, {}).get("state") == "failed"]
    return 1 if failed else 0


def status():
    st = load_state()
    print(f"{'job':26s} {'state':12s} {'validation':22s} reasons")
    for jid in JOBS:
        ok, s, why = validated(jid)
        print(f"{jid:26s} {st.get(jid, {}).get('state', '-'):12s} {('OK ' if ok else '') + str(s):22s} "
              f"{'; '.join(why)[:110]}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run")
    r.add_argument("--max-jobs", type=int, default=2)
    r.add_argument("--only", default=None)
    r.add_argument("--retry-failed", action="store_true")
    sp.add_parser("status")
    s = sp.add_parser("show")
    s.add_argument("job")
    a = p.parse_args(argv)
    if a.cmd == "status":
        status()
        return 0
    if a.cmd == "show":
        j = JOBS[a.job]
        print(json.dumps(dict(command=command(a.job) if j.get("kind") != "legacy_audit" else "legacy audit",
                              outdir=outdir(a.job), deps=j["deps"], threads=j["threads"], mem_mb=j["mem_mb"],
                              allow_provisional=j.get("allow_provisional")), indent=1))
        return 0
    return run(a.max_jobs, a.only, retry_failed=a.retry_failed)


if __name__ == "__main__":
    sys.exit(main())
