"""Resource budget: cgroup-aware limits, admission, stale reservations, disk and thread checks, and the
budgeted runner. Uses small deterministic fixtures (fake /proc and cgroup trees); nothing is pushed to
the real container limits."""

import importlib
import json
import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))


def _fixture(tmp_path, limit_mb=10000, rss_mb=3000, cache_mb=1000, inactive_mb=800, ver=1):
    proc, cg = tmp_path / "proc", tmp_path / "cg"
    (proc / "self").mkdir(parents=True)
    (proc / "meminfo").write_text("MemTotal:       16000000 kB\n")
    mb = 1024 * 1024
    if ver == 1:
        (proc / "self" / "cgroup").write_text("4:memory:/a/b\n1:cpu:/\n")
        d = cg / "memory" / "a" / "b"
        d.mkdir(parents=True)
        (cg / "memory" / "memory.limit_in_bytes").write_text(str(1 << 62))
        (cg / "memory" / "a" / "memory.limit_in_bytes").write_text(str(1 << 62))
        (d / "memory.limit_in_bytes").write_text(str(limit_mb * mb))
        (d / "memory.usage_in_bytes").write_text(str((rss_mb + cache_mb) * mb))
        (d / "memory.stat").write_text(f"total_rss {rss_mb * mb}\ntotal_cache {cache_mb * mb}\n"
                                       f"total_inactive_file {inactive_mb * mb}\ntotal_mapped_file {200 * mb}\n")
    else:
        (proc / "self" / "cgroup").write_text("0::/x\n")
        d = cg / "x"
        d.mkdir(parents=True)
        (d / "memory.max").write_text(str(limit_mb * mb))
        (d / "memory.current").write_text(str((rss_mb + cache_mb) * mb))
        (d / "memory.stat").write_text(f"anon {rss_mb * mb}\nfile {cache_mb * mb}\ninactive_file {inactive_mb * mb}\n"
                                       f"file_mapped {200 * mb}\n")
    return proc, cg


def _load(monkeypatch, tmp_path, proc, cg):
    monkeypatch.setenv("ALGAE_PROC", str(proc))
    monkeypatch.setenv("ALGAE_CGROUP_FS", str(cg))
    monkeypatch.setenv("ALGAE_BUDGET_LEDGER", str(tmp_path / "ledger" / "budget.json"))
    import resources
    return importlib.reload(resources)


@pytest.mark.parametrize("ver", [1, 2])
def test_cgroup_limit_not_host_memory(tmp_path, monkeypatch, ver):
    proc, cg = _fixture(tmp_path, limit_mb=10000, ver=ver)
    R = _load(monkeypatch, tmp_path, proc, cg)
    ms = R.memory_status()
    assert abs(ms["limit_mb"] - 10000) < 1 and ms["host_total_mb"] > 15000
    assert abs(ms["unreclaimable_mb"] - (3000 + 200)) < 1           # anon rss + mapped file


def test_admission_counts_reservations_and_refuses_overcommit(tmp_path, monkeypatch):
    proc, cg = _fixture(tmp_path, limit_mb=10000, rss_mb=3000)
    R = _load(monkeypatch, tmp_path, proc, cg)
    monkeypatch.setattr(R, "tree_rss_mb", lambda pid: 0.0)            # reserved jobs not yet grown
    ok, why, _ = R.try_reserve("chem", os.getpid(), 4000, threads=2, cores=4)
    assert ok
    ok, why, _ = R.try_reserve("analysis", os.getpid(), 2500, threads=1, cores=4)
    assert not ok and any("memory" in w for w in why)                # 3200 + 4000 + 2500 > 10000 - 1200
    ok, _, _ = R.try_reserve("small", os.getpid(), 1000, threads=1, cores=4)
    assert ok
    ok, why, _ = R.try_reserve("t", os.getpid(), 10, threads=2, cores=4)
    assert not ok and any("threads" in w for w in why)


def test_stale_reservations_are_pruned(tmp_path, monkeypatch):
    proc, cg = _fixture(tmp_path, limit_mb=10000, rss_mb=1000)
    R = _load(monkeypatch, tmp_path, proc, cg)
    os.makedirs(os.path.dirname(R.LEDGER), exist_ok=True)
    json.dump({"dead": dict(pid=999999999, mem_mb=8000, threads=4)}, open(R.LEDGER, "w"))
    monkeypatch.setattr(R, "tree_rss_mb", lambda pid: 0.0)
    ok, why, _ = R.try_reserve("new", os.getpid(), 4000, threads=2, cores=4)
    assert ok and "dead" not in R.snapshot()["reservations"]


def test_disk_admission(tmp_path, monkeypatch):
    proc, cg = _fixture(tmp_path)
    R = _load(monkeypatch, tmp_path, proc, cg)
    free = R.shutil.disk_usage(str(tmp_path)).free / 1e9
    ok, why, _ = R.try_reserve("big", os.getpid(), 10, scratch_gb=free + 1, scratch_dir=str(tmp_path), cores=64)
    assert not ok and any("disk" in w for w in why)


def test_thread_env_caps_all_pools():
    import resources as R
    env = R.thread_env(2, base={})
    assert all(env[k] == "2" for k in R.THREAD_VARS)


def test_run_budgeted_refuses_and_runs(tmp_path, monkeypatch):
    proc, cg = _fixture(tmp_path, limit_mb=10000, rss_mb=3000)
    env = dict(os.environ, ALGAE_PROC=str(proc), ALGAE_CGROUP_FS=str(cg),
               ALGAE_BUDGET_LEDGER=str(tmp_path / "ledger" / "budget.json"))
    rb = os.path.join(HERE, "..", "run_budgeted.py")
    out = tmp_path / "o.txt"
    # too big -> not admitted, command not run, exit 75
    r = subprocess.run([sys.executable, rb, "--name", "x", "--mem-mb", "9000", "--scratch-root", str(tmp_path / "s"),
                        "--", sys.executable, "-c", f"open({str(out)!r},'w').write('ran')"], env=env)
    assert r.returncode == 75 and not out.exists()
    # fits -> runs with thread caps and a private scratch, reservation released afterwards
    code = (f"import os; open({str(out)!r},'w').write(os.environ['OMP_NUM_THREADS']+'|'+os.environ['TMPDIR'])")
    r = subprocess.run([sys.executable, rb, "--name", "y", "--mem-mb", "500", "--threads", "1",
                        "--scratch-root", str(tmp_path / "s"), "--", sys.executable, "-c", code], env=env)
    assert r.returncode == 0 and out.read_text().startswith("1|") and out.read_text().endswith("/y")
    assert json.load(open(tmp_path / "ledger" / "budget.json")) == {}
    # the RLIMIT_AS backstop turns a runaway allocation into a failure of that process only
    r = subprocess.run([sys.executable, rb, "--name", "z", "--mem-mb", "200", "--as-factor", "2",
                        "--scratch-root", str(tmp_path / "s"), "--", sys.executable, "-c",
                        "b = bytearray(2 * 1024**3)"], env=env, capture_output=True)
    assert r.returncode != 0 and b"MemoryError" in r.stderr


def test_background_class_beyond_cores_but_memory_strict(tmp_path, monkeypatch):
    proc, cg = _fixture(tmp_path, limit_mb=10000, rss_mb=3000)
    R = _load(monkeypatch, tmp_path, proc, cg)
    monkeypatch.setattr(R, "tree_rss_mb", lambda pid: 0.0)
    assert R.try_reserve("chem", os.getpid(), 3000, threads=4, cores=4)[0]          # all cores reserved
    ok, why, _ = R.try_reserve("fg", os.getpid(), 100, threads=1, cores=4)
    assert not ok and any("threads" in w for w in why)
    assert R.try_reserve("bg1", os.getpid(), 500, threads=1, cores=4, background=True)[0]
    assert R.try_reserve("bg2", os.getpid(), 100, threads=1, cores=4, background=True)[0]
    ok, why, _ = R.try_reserve("bg2b", os.getpid(), 100, threads=1, cores=4, background=True)
    assert not ok and any("background" in w for w in why)                            # BACKGROUND_SLOTS = 2
    R.release("bg1")
    R.release("bg2")
    ok, why, _ = R.try_reserve("bg3", os.getpid(), 4000, threads=1, cores=4, background=True)
    assert not ok and any("memory" in w for w in why)                                # memory unchanged


def test_run_budgeted_background_runs_niced(tmp_path):
    proc, cg = _fixture(tmp_path, limit_mb=10000, rss_mb=3000)
    env = dict(os.environ, ALGAE_PROC=str(proc), ALGAE_CGROUP_FS=str(cg),
               ALGAE_BUDGET_LEDGER=str(tmp_path / "ledger" / "budget.json"))
    rb = os.path.join(HERE, "..", "run_budgeted.py")
    out = tmp_path / "n.txt"
    r = subprocess.run([sys.executable, rb, "--name", "b", "--mem-mb", "300", "--background",
                        "--scratch-root", str(tmp_path / "s"), "--", sys.executable, "-c",
                        f"import os; open({str(out)!r},'w').write(str(os.nice(0)))"], env=env)
    assert r.returncode == 0 and out.read_text() == "19"


def test_short_class_limits(tmp_path, monkeypatch):
    proc, cg = _fixture(tmp_path, limit_mb=10000, rss_mb=3000)
    R = _load(monkeypatch, tmp_path, proc, cg)
    monkeypatch.setattr(R, "tree_rss_mb", lambda pid: 0.0)
    assert R.try_reserve("chem", os.getpid(), 3000, threads=4, cores=4)[0]
    assert not R.try_reserve("big", os.getpid(), 900, threads=1, cores=4, short=True)[0]
    assert R.try_reserve("s1", os.getpid(), 300, threads=1, cores=4, short=True)[0]
    ok, why, _ = R.try_reserve("s2", os.getpid(), 300, threads=1, cores=4, short=True)
    assert not ok and any("short" in w for w in why)


def test_run_budgeted_short_is_killed_at_limit(tmp_path, monkeypatch):
    proc, cg = _fixture(tmp_path, limit_mb=10000, rss_mb=3000)
    env = dict(os.environ, ALGAE_PROC=str(proc), ALGAE_CGROUP_FS=str(cg),
               ALGAE_BUDGET_LEDGER=str(tmp_path / "ledger" / "budget.json"))
    rb = os.path.join(HERE, "..", "run_budgeted.py")
    code = "import resources as R; R.SHORT_MAX_S = 1; import sys; sys.argv = sys.argv[1:]; import run_budgeted as B; sys.exit(B.main(sys.argv))"
    r = subprocess.run([sys.executable, "-c", code, "--name", "s", "--mem-mb", "100", "--short", "--scratch-root",
                        str(tmp_path / "s"), "--", sys.executable, "-c", "import time; time.sleep(30)"],
                       env=dict(env, PYTHONPATH=os.path.join(HERE, "..")), timeout=60)
    assert r.returncode == 124


def _idproc(tmp_path, pid, start_ticks, boot="boot-A", cmd="python3 run_phase1.py --outdir /x/JOB", btime=1000):
    d = tmp_path / "idproc"
    (d / "sys" / "kernel" / "random").mkdir(parents=True, exist_ok=True)
    (d / "sys" / "kernel" / "random" / "boot_id").write_text(boot + "\n")
    (d / "stat").write_text(f"cpu 1 2 3\nbtime {btime}\n")
    (d / str(pid)).mkdir(exist_ok=True)
    fields = ["S"] + ["0"] * 18 + [str(start_ticks)] + ["0"] * 10
    (d / str(pid) / "stat").write_text(f"{pid} (python3) " + " ".join(fields) + "\n")
    (d / str(pid) / "cmdline").write_bytes(cmd.replace(" ", "\0").encode())
    return d


def test_identity_rejects_reboot_and_pid_reuse(tmp_path, monkeypatch):
    pid = os.getpid()                                  # alive, so only identity can reject it
    d = _idproc(tmp_path, pid, start_ticks=500)
    monkeypatch.setenv("ALGAE_PROC_IDENTITY", str(d))
    import resources
    R = importlib.reload(resources)
    ident = R.process_identity(pid)
    assert ident["boot_id"] == "boot-A" and ident["start_epoch"] == 1000 + 500 / os.sysconf("SC_CLK_TCK")
    assert R.identity_ok(ident)
    # reboot: same pid alive, different boot id
    _idproc(tmp_path, pid, start_ticks=500, boot="boot-B")
    assert not R.identity_ok(ident)
    # pid reuse in the same boot: same pid, later start time
    _idproc(tmp_path, pid, start_ticks=90000)
    assert not R.identity_ok(ident)
    # same start time, different command line
    _idproc(tmp_path, pid, start_ticks=500, cmd="python3 other.py")
    assert not R.identity_ok(ident)
    # legacy record (no identity): accepted only if the process started before the recorded launch
    _idproc(tmp_path, pid, start_ticks=500)
    assert R.identity_ok(dict(pid=pid), started=1000 + 600, must_contain="/x/JOB")
    assert not R.identity_ok(dict(pid=pid), started=1000 - 100, must_contain="/x/JOB")   # process younger than the record -> reused
    assert not R.identity_ok(dict(pid=pid), started=1000 + 600, must_contain="/x/OTHER")


def test_ledger_prunes_reservation_of_reused_pid(tmp_path, monkeypatch):
    proc, cg = _fixture(tmp_path, limit_mb=10000, rss_mb=3000)
    pid = os.getpid()
    d = _idproc(tmp_path, pid, start_ticks=500)
    monkeypatch.setenv("ALGAE_PROC_IDENTITY", str(d))
    R = _load(monkeypatch, tmp_path, proc, cg)
    monkeypatch.setattr(R, "tree_rss_mb", lambda p: 0.0)
    assert R.try_reserve("chem", pid, 4000, threads=2, cores=4)[0]
    assert "chem" in R.snapshot()["reservations"]
    _idproc(tmp_path, pid, start_ticks=500, boot="boot-B")          # container rebooted, pid reused
    assert "chem" not in R.snapshot()["reservations"]
