"""Job runner: completion requires a validated marker; failures propagate; dependencies block;
interrupted jobs are relaunched."""

import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import jobs  # noqa: E402

GOOD = ("import completion,sys,os,json; d=sys.argv[1]; os.makedirs(d,exist_ok=True); "
        "p=os.path.join(d,'summary.json'); json.dump({'tddft':{}},open(p,'w')); "
        "completion.write_marker(d,{'tddft':{}},[p],True,{'status':'converged'})")


@pytest.fixture
def reg(tmp_path, monkeypatch):
    monkeypatch.setattr(jobs, "V2", str(tmp_path))
    monkeypatch.setattr(jobs, "STATE", str(tmp_path / "state.json"))
    monkeypatch.setattr(jobs, "external_phase1_threads", lambda own: 0)    # isolate from other processes
    monkeypatch.setattr(jobs, "mem_available_mb", lambda: 1e6)
    monkeypatch.setattr(jobs, "mem_total_mb", lambda: 1e6)
    monkeypatch.setattr(jobs, "external_phase1_rss_mb", lambda own: 0)
    monkeypatch.setattr(jobs, "admit", lambda jid, j: (True, []))       # isolate from the real shared budget
    monkeypatch.setattr(jobs, "register_running", lambda st: None)
    monkeypatch.setattr(jobs.RS, "update_pid", lambda name, pid: None)
    monkeypatch.setattr(jobs.RS, "release", lambda name: None)
    py = sys.executable

    def job(code, deps=()):
        return lambda d: dict(cmd=[py, "-c", code, d], deps=list(deps), threads=1, mem_mb=10)
    defs = {
        "ok": job(GOOD),
        "fails": job("import sys; sys.exit(1)"),
        "summary_only": job("import sys,os,json; d=sys.argv[1]; os.makedirs(d,exist_ok=True); "
                            "json.dump({'x':1},open(os.path.join(d,'summary.json'),'w'))"),
        "after_ok": job(GOOD, deps=["ok"]),
        "after_fail": job(GOOD, deps=["fails"]),
    }
    reg = {k: f(str(tmp_path / k)) for k, f in defs.items()}
    monkeypatch.setattr(jobs, "JOBS", reg)
    return tmp_path


def test_queue_failure_and_dependencies(reg):
    rc = jobs.run(max_jobs=2, poll=0.2)
    st = json.load(open(reg / "state.json"))
    assert rc == 1                                            # a failure is never reported as success
    assert st["ok"]["state"] == "complete" and st["after_ok"]["state"] == "complete"
    assert st["fails"]["state"] == "failed"
    assert st["summary_only"]["state"] == "failed"            # nonempty summary.json is not completion
    assert "after_fail" not in st or st["after_fail"].get("state") != "complete"


def test_interrupted_job_is_relaunched(reg):
    jobs.save_state({"ok": dict(state="running", pid=999999)})
    rc = jobs.run(max_jobs=1, only="ok", poll=0.2)
    assert rc == 0 and json.load(open(reg / "state.json"))["ok"]["state"] == "complete"


def test_scratch_is_per_job_and_cleaned_before_relaunch(tmp_path, monkeypatch):
    import jobs
    monkeypatch.setattr(jobs, "outdir", lambda jid: str(tmp_path / jid))
    d = jobs.clean_scratch("X")
    stale = os.path.join(d, "orphan_cderi.h5")
    open(stale, "w").write("x" * 100)
    assert os.path.isfile(stale)
    d2 = jobs.clean_scratch("X")                       # relaunch: orphaned DF tensors of a killed run go
    assert d2 == d and os.listdir(d2) == []
    monkeypatch.setitem(jobs.JOBS, "Y", dict(args=["--level", "level2", "--skip-opt"], threads=1, mem_mb=1000, deps=[]))
    cmd = jobs.command("Y")
    assert cmd[cmd.index("--td-chunk") + 1] == str(jobs.TD_CHUNK)


def test_runner_waits_when_budget_refuses(tmp_path, monkeypatch):
    import jobs
    monkeypatch.setattr(jobs, "V2", str(tmp_path))
    monkeypatch.setattr(jobs, "STATE", str(tmp_path / "state.json"))
    monkeypatch.setattr(jobs, "JOBS", {"A": dict(cmd=[sys.executable, "-c", "pass"], deps=[], threads=1, mem_mb=100)})
    monkeypatch.setattr(jobs, "outdir", lambda jid: str(tmp_path / jid))
    monkeypatch.setattr(jobs, "validated", lambda jid: (False, "missing", ["no marker"]))
    monkeypatch.setattr(jobs, "register_running", lambda st: None)
    calls = []

    def refuse(jid, j):
        calls.append(jid)
        if len(calls) > 2:                     # stop the loop: pretend the job was tried
            raise KeyboardInterrupt
        return False, ["memory: projected over limit"]
    monkeypatch.setattr(jobs, "admit", refuse)
    with pytest.raises(KeyboardInterrupt):
        jobs.run(max_jobs=1, poll=0)
    assert calls == ["A", "A", "A"] and not os.path.exists(tmp_path / "A" / "job.log")   # never launched


def test_stale_running_record_with_reused_pid_is_relaunched(reg):
    """After a reboot the recorded pid can belong to an unrelated live process: it must not be adopted."""
    import time
    me = os.getpid()                              # alive, but not the job (cmdline lacks the job's outdir)
    jobs.save_state({"ok": dict(state="running", pid=me, started=time.time() - 1e5)})
    rc = jobs.run(max_jobs=1, only="ok", poll=0.2)
    assert rc == 0 and json.load(open(reg / "state.json"))["ok"]["state"] == "complete"


def test_identity_mismatch_marks_interrupted(reg, monkeypatch):
    ident = jobs.RS.process_identity(os.getpid())
    ident = dict(ident, boot_id="a-previous-boot")
    rec = dict(state="running", pid=os.getpid(), started=0.0, ident=ident)
    assert not jobs.job_alive("ok", rec)
    rec["ident"] = jobs.RS.process_identity(os.getpid())
    monkeypatch.setattr(jobs, "outdir", lambda jid: "")        # identity matches; outdir check trivially true
    assert jobs.job_alive("ok", rec)


def test_converged_optimisation_is_not_repeated(tmp_path, monkeypatch):
    od = tmp_path / "OPTJOB"
    od.mkdir()
    (od / "opt_B3LYP_final.xyz").write_text("1\nH\nH 0 0 0\n")
    (od / "opt_B3LYP_opt_record.json").write_text(json.dumps(dict(converged=True, optimizer="geomeTRIC",
                                                                  thresholds={"convergence_grms": 3e-4})))
    monkeypatch.setattr(jobs, "outdir", lambda jid: str(od))
    monkeypatch.setattr(jobs, "JOBS", {"OPTJOB": dict(args=["--level", "level1", "--nstates", "30"], deps=[],
                                                      threads=1, mem_mb=10)})
    cmd = jobs.command("OPTJOB")
    assert "--skip-opt" in cmd and cmd[cmd.index("--start-xyz") + 1] == str(od / "opt_B3LYP_final.xyz")
    side = json.load(open(od / "opt_B3LYP_final.xyz.geometry.json"))
    assert side["status"] == "converged" and len(side["xyz_sha256"]) == 64
    # an unconverged optimisation keeps resuming from its last geometry
    (od / "opt_B3LYP_opt_record.json").write_text(json.dumps(dict(converged=False)))
    (od / "opt_B3LYP_last.xyz").write_text("1\nH\nH 0 0 0.1\n")
    cmd = jobs.command("OPTJOB")
    assert "--skip-opt" not in cmd and cmd[cmd.index("--start-xyz") + 1] == str(od / "opt_B3LYP_last.xyz")
