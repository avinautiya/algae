#!/usr/bin/env python
"""
Run a command only after the shared resource budget admits it (common/resources.py).

    python common/run_budgeted.py --name heldout --mem-mb 2000 --threads 1 --scratch-gb 1 \
        [--wait 3600] -- python3 phase4/heldout.py ...

  * waits (polling) up to --wait seconds for admission, else exits 75 (EX_TEMPFAIL) without running;
  * sets OMP/OpenBLAS/MKL/NumExpr/VecLib/BLIS thread counts to --threads;
  * --background: admitted beyond the core count (at most resources.BACKGROUND_SLOTS threads in total) and
    run at nice 19, so it only uses CPU the reserved jobs leave idle; memory admission is unchanged;
  * applies RLIMIT_AS = --as-factor x mem (hard backstop: the process fails with MemoryError instead of
    pushing the container's cgroup into an OOM kill of other jobs);
  * per-run scratch TMPDIR under --scratch-root/<name>, removed after a successful run, kept on failure
    for inspection; the reservation is released when the process exits.
"""

from __future__ import annotations

import argparse
import os
import resource
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import resources as RS  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--name", required=True)
    p.add_argument("--mem-mb", type=float, required=True)
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--scratch-gb", type=float, default=0.5)
    p.add_argument("--scratch-root", default=os.path.join(RS.ROOT, "phase1", "results", "resources", "scratch"))
    p.add_argument("--wait", type=float, default=0.0, help="seconds to wait for admission")
    p.add_argument("--as-factor", type=float, default=2.5, help="RLIMIT_AS = factor x mem-mb (virtual > resident)")
    p.add_argument("--background", action="store_true",
                   help=f"background class: beyond the core count (max {RS.BACKGROUND_SLOTS} thread), nice {RS.BACKGROUND_NICE}")
    p.add_argument("--short", action="store_true",
                   help=f"short class: <= {RS.SHORT_MAX_MB} MB, 1 thread, killed after {RS.SHORT_MAX_S} s, one at a time")
    p.add_argument("cmd", nargs=argparse.REMAINDER)
    a = p.parse_args(argv)
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        p.error("no command")
    deadline = time.time() + a.wait
    while True:
        ok, reasons, snap = RS.try_reserve(a.name, os.getpid(), a.mem_mb, a.threads, a.scratch_gb, background=a.background,
                                         short=a.short)
        if ok:
            break
        if time.time() >= deadline:
            print(f"[budget] {a.name} NOT admitted: {'; '.join(reasons)}", file=sys.stderr, flush=True)
            return 75
        time.sleep(30)
    scratch = os.path.join(a.scratch_root, a.name)
    shutil.rmtree(scratch, ignore_errors=True)
    os.makedirs(scratch, exist_ok=True)
    env = RS.thread_env(a.threads)
    env.update(TMPDIR=scratch, PYSCF_TMPDIR=scratch)
    limit = int(a.as_factor * a.mem_mb * 1024 * 1024)

    def pre():
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
        if a.background or a.short:
            os.nice(RS.BACKGROUND_NICE - os.nice(0))
    print(f"[budget] {a.name} admitted: {a.mem_mb:.0f} MB, {a.threads} threads; projected "
          f"{snap['projected_mb']:.0f}/{snap['memory']['limit_mb']:.0f} MB", file=sys.stderr, flush=True)
    try:
        proc = subprocess.Popen(cmd, env=env, preexec_fn=pre)
        RS.update_pid(a.name, proc.pid)
        try:
            rc = proc.wait(timeout=RS.SHORT_MAX_S if a.short else None)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            print(f"[budget] {a.name} exceeded the short-class limit of {RS.SHORT_MAX_S} s: killed", file=sys.stderr)
            rc = 124
    finally:
        RS.release(a.name)
    if rc == 0:
        shutil.rmtree(scratch, ignore_errors=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
