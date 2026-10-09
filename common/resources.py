"""
Shared compute-resource budget for every expensive process in this repository (chemistry jobs,
calibration, emulator builds, scoring).

Why: two OOM kills and a disk-full failure (docs/repair_ledger.md P1-ENV-2/3) happened because
admission used the HOST memory (16 GB) and a fixed reserve, while the governing limit is the memory
cgroup of this container (14.35 GB, `claude-code-bash`), and analysis processes were not budgeted at all.

Model:
  * memory: the effective limit is the minimum `memory.limit_in_bytes` (v1) / `memory.max` (v2) along the
    process's cgroup path. Usage that cannot be reclaimed = rss + (cache - inactive_file).
  * reservations: every budgeted process registers {mem_mb, threads, scratch_gb, pid} in a shared,
    file-locked ledger. A live reservation counts at max(reserved, its process tree's current RSS);
    unregistered processes count at their measured RSS (part of cgroup usage).
  * admission: (unreclaimable usage not covered by reservations) + sum(live reservations) + request
    <= limit - margin;  threads <= cores;  disk free >= scratch + disk margin.
  * thread pools: OMP/OpenBLAS/MKL/NumExpr/VecLib/BLIS thread counts are set for the admitted process.

`nice` is not a memory limit; RLIMIT_AS is applied as a hard backstop in run_budgeted.
"""

from __future__ import annotations

import contextlib
import fcntl
import json
import os
import shutil
import time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LEDGER = os.environ.get("ALGAE_BUDGET_LEDGER", os.path.join(ROOT, "phase1", "results", "resources", "budget.json"))
CGROUP_FS = os.environ.get("ALGAE_CGROUP_FS", "/sys/fs/cgroup")
PROC = os.environ.get("ALGAE_PROC", "/proc")
MEM_MARGIN_MB = 1200            # kernel, page tables, transient spikes
DISK_MARGIN_GB = 5.0
# Background class: analysis admitted beyond the core count, at most BACKGROUND_SLOTS threads in total,
# run at nice 19 (CFS weight 15 vs 1024: < 2 % of a contended core). It uses idle CPU only; memory
# admission is unchanged (memory, not CPU, is what can kill chemistry).
BACKGROUND_SLOTS = 1
BACKGROUND_NICE = 19
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")
_UNLIMITED = 1 << 60


# --------------------------------------------------------------------------- memory
def _cgroup_paths():
    """(version, relative path) of this process's memory cgroup."""
    with open(os.path.join(PROC, "self", "cgroup")) as fh:
        lines = fh.read().splitlines()
    for ln in lines:
        _, ctrl, path = ln.split(":", 2)
        if "memory" in ctrl.split(","):
            return 1, path
    for ln in lines:
        if ln.startswith("0::"):
            return 2, ln.split("::", 1)[1]
    return None, None


def _read_int(p):
    try:
        v = open(p).read().strip()
    except OSError:
        return None
    if v == "max":
        return _UNLIMITED
    try:
        return int(v)
    except ValueError:
        return None


def memory_status() -> dict:
    """Effective cgroup memory limit and usage in MB (host values if no cgroup is found)."""
    ver, path = _cgroup_paths()
    limit, usage, stat = None, None, {}
    if ver == 1:
        base = os.path.join(CGROUP_FS, "memory")
        parts = [p for p in path.strip("/").split("/") if p]
        limits = []
        for k in range(len(parts), -1, -1):
            d = os.path.join(base, *parts[:k])
            v = _read_int(os.path.join(d, "memory.limit_in_bytes"))
            if v is not None:
                limits.append(v)
        here = os.path.join(base, *parts)
        limit = min(limits) if limits else None
        usage = _read_int(os.path.join(here, "memory.usage_in_bytes"))
        stat_file = os.path.join(here, "memory.stat")
    elif ver == 2:
        parts = [p for p in path.strip("/").split("/") if p]
        limits = [v for k in range(len(parts), -1, -1)
                  if (v := _read_int(os.path.join(CGROUP_FS, *parts[:k], "memory.max"))) is not None]
        here = os.path.join(CGROUP_FS, *parts)
        limit = min(limits) if limits else None
        usage = _read_int(os.path.join(here, "memory.current"))
        stat_file = os.path.join(here, "memory.stat")
    else:
        stat_file = None
    if stat_file and os.path.isfile(stat_file):
        for ln in open(stat_file):
            k, v = ln.split()
            stat[k] = int(v)
    host_total = None
    for ln in open(os.path.join(PROC, "meminfo")):
        if ln.startswith("MemTotal:"):
            host_total = int(ln.split()[1]) * 1024
    if limit is None or limit >= _UNLIMITED // 2:
        limit = host_total
    rss = stat.get("total_rss", stat.get("rss", stat.get("anon")))
    mapped = stat.get("total_mapped_file", stat.get("mapped_file", stat.get("file_mapped", 0)))
    shmem = stat.get("total_shmem", stat.get("shmem", 0))
    if rss is None:
        rss = usage or 0
    # memory the kernel cannot drop under pressure: anonymous + mapped file + shmem. Unmapped page cache
    # (active or inactive) is reclaimed before an OOM kill, so it is not counted (the margin covers
    # kernel-side usage and transient spikes)
    unreclaimable = rss + mapped + shmem
    mb = 1024 * 1024
    return dict(cgroup_version=ver, cgroup_path=path, limit_mb=limit / mb, usage_mb=(usage or 0) / mb,
                rss_mb=rss / mb, unreclaimable_mb=unreclaimable / mb, host_total_mb=host_total / mb)


def tree_rss_mb(pid: int) -> float:
    """RSS of a process and its descendants (MB); 0 if gone."""
    total, todo, seen = 0, [int(pid)], set()
    children = {}
    for d in os.listdir(PROC):
        if d.isdigit():
            try:
                pp = int(open(os.path.join(PROC, d, "stat")).read().rsplit(")", 1)[1].split()[1])
                children.setdefault(pp, []).append(int(d))
            except (OSError, IndexError, ValueError):
                pass
    while todo:
        p = todo.pop()
        if p in seen:
            continue
        seen.add(p)
        try:
            for ln in open(os.path.join(PROC, str(p), "status")):
                if ln.startswith("VmRSS:"):
                    total += int(ln.split()[1])
        except OSError:
            continue
        todo += children.get(p, [])
    return total / 1024.0


def pid_alive(pid) -> bool:
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, TypeError, ValueError):
        return False


# --------------------------------------------------------------------------- ledger
@contextlib.contextmanager
def _locked():
    os.makedirs(os.path.dirname(LEDGER), exist_ok=True)
    with open(LEDGER + ".lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            try:
                led = json.load(open(LEDGER))
            except (OSError, ValueError):
                led = {}
            # prune reservations of dead processes
            led = {k: v for k, v in led.items() if pid_alive(v.get("pid"))}
            yield led
            tmp = LEDGER + f".tmp{os.getpid()}"
            with open(tmp, "w") as fh:
                json.dump(led, fh, indent=1)
            os.replace(tmp, LEDGER)
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


def _assess(led, mem_mb, threads, scratch_gb, scratch_dir, cores=None, background=False):
    ms = memory_status()
    reserved = 0.0
    covered = 0.0
    for v in led.values():
        cur = tree_rss_mb(v["pid"])
        reserved += max(v["mem_mb"], cur)
        covered += cur
    other = max(ms["unreclaimable_mb"] - covered, 0.0)        # unregistered processes + kernel-side usage
    projected = other + reserved + mem_mb
    cores = cores or os.cpu_count() or 1
    threads_used = sum(v.get("threads", 1) for v in led.values() if not v.get("background"))
    bg_used = sum(v.get("threads", 1) for v in led.values() if v.get("background"))
    free_gb = shutil.disk_usage(scratch_dir or ROOT).free / 1e9
    reasons = []
    if projected > ms["limit_mb"] - MEM_MARGIN_MB:
        reasons.append(f"memory: projected {projected:.0f} MB > limit {ms['limit_mb']:.0f} - margin {MEM_MARGIN_MB}")
    if background:
        if bg_used + threads > BACKGROUND_SLOTS:
            reasons.append(f"background threads: {bg_used} reserved + {threads} > {BACKGROUND_SLOTS} slot(s)")
    elif threads_used + threads > cores:
        reasons.append(f"threads: {threads_used} reserved + {threads} > {cores} cores")
    if free_gb < scratch_gb + DISK_MARGIN_GB:
        reasons.append(f"disk: {free_gb:.1f} GB free < {scratch_gb} + {DISK_MARGIN_GB} GB margin")
    return reasons, dict(memory=ms, reserved_mb=reserved, unregistered_mb=other, projected_mb=projected,
                         threads_reserved=threads_used, disk_free_gb=free_gb)


def try_reserve(name, pid, mem_mb, threads=1, scratch_gb=0.0, scratch_dir=None, cores=None, background=False):
    """Atomically admit and register a reservation. Returns (admitted, reasons, snapshot).
    background=True: admitted beyond the core count within BACKGROUND_SLOTS; the caller MUST run it at
    nice BACKGROUND_NICE (run_budgeted does)."""
    with _locked() as led:
        reasons, snap = _assess(led, mem_mb, threads, scratch_gb, scratch_dir, cores, background)
        if not reasons:
            led[name] = dict(pid=int(pid), mem_mb=float(mem_mb), threads=int(threads), scratch_gb=float(scratch_gb),
                             started=time.time(), background=bool(background))
        return not reasons, reasons, snap


def update_pid(name, pid):
    with _locked() as led:
        if name in led:
            led[name]["pid"] = int(pid)


def release(name):
    with _locked() as led:
        led.pop(name, None)


def snapshot():
    with _locked() as led:
        return dict(reservations=led, memory=memory_status())


def thread_env(threads: int, base=None) -> dict:
    env = dict(os.environ if base is None else base)
    for k in THREAD_VARS:
        env[k] = str(int(threads))
    return env


# --------------------------------------------------------------------------- estimates
# Conservative peak-memory estimates (MB) from observed runs on this container (P1-ENV-2/3, Phase 4 runs):
ESTIMATES_MB = {
    "phase1_level2_td": 7500,       # observed RSS 6.4-7.1 GB
    "phase1_level2_opt": 6000,      # observed 5.5 GB
    "phase1_level1": 3000,          # observed 1.9 GB
    "phase1_fe": 4500,              # observed 3.1 GB (still growing at kill)
    "tddft_calibration_ar1": 1500,  # stage-2 grids (x1)
    "emulator_build_2workers": 1500,
    "heldout_scoring": 2000,
    "surrogate_benchmark": 1500,
    "scene_product": 2500,
}
