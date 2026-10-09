#!/usr/bin/env python
"""Portable execution / recovery manifest for the Phase 1 chemistry queue (run on a stable machine).

Writes records/portable_chemistry/manifest.json: for every incomplete job the exact command, input files and
checkpoints with sha256, declared memory/threads/scratch, code hashes and pinned dependencies. Nothing is
launched. See docs/portable_chemistry_package.md for the procedure.

    python3 phase1/portable_package.py
"""
import glob
import hashlib
import json
import os
import platform
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
import jobs  # noqa: E402


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rel(p):
    return os.path.relpath(os.path.abspath(p), os.path.abspath(ROOT))


def main():
    out = {}
    for jid, j in jobs.JOBS.items():
        if j.get("kind") == "legacy_audit" or jobs.validated(jid)[0]:
            continue
        cmd = jobs.command(jid) if j.get("args") is not None or j.get("geometry_from") else None
        inputs = {}
        if cmd:
            for flag in ("--start-xyz", "--import-guess"):
                if flag in cmd and os.path.isfile(cmd[cmd.index(flag) + 1]):
                    p = cmd[cmd.index(flag) + 1]
                    inputs[rel(p)] = sha(p)
                    if os.path.isfile(p + ".geometry.json"):
                        inputs[rel(p + ".geometry.json")] = sha(p + ".geometry.json")
        od = jobs.outdir(jid)
        ckpts = {rel(p): dict(sha256=sha(p), bytes=os.path.getsize(p))
                 for p in sorted(glob.glob(os.path.join(od, "td_ckpt", "*", "*")) + glob.glob(os.path.join(od, "*.chk"))
                                 + glob.glob(os.path.join(od, "opt_*")))
                 if os.path.isfile(p) and not p.endswith(".lock")}
        out[jid] = dict(command=[rel(c) if os.path.exists(c) else c for c in cmd] if cmd else None,
                        deps=j["deps"], threads=j["threads"], declared_peak_mem_mb=j["mem_mb"],
                        scratch_gb=8.0 if j["mem_mb"] >= 5000 else 2.0, inputs=inputs, resume_files=ckpts,
                        exploratory=j.get("exploratory"), allow_provisional=j.get("allow_provisional"))
    code = {rel(p): sha(p) for p in sorted(glob.glob(os.path.join(HERE, "*.py")))}
    man = dict(jobs=out, code=code, requirements=open(os.path.join(ROOT, "requirements.txt")).read(),
               python=platform.python_version(), launch_order=list(jobs.JOBS),
               note="Copy the repository including phase1/results/v2/<job>/ (td_ckpt, *.chk, opt_*) and "
                    "phase1/results/level2/; verify sha256; run `python3 phase1/jobs.py run --max-jobs N`.")
    d = os.path.join(ROOT, "records", "portable_chemistry")
    os.makedirs(d, exist_ok=True)
    json.dump(man, open(os.path.join(d, "manifest.json"), "w"), indent=1)
    print(json.dumps({k: dict(mem=v["declared_peak_mem_mb"], n_resume_files=len(v["resume_files"])) for k, v in out.items()}, indent=1))


if __name__ == "__main__":
    main()
