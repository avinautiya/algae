"""Campaign inventory from case lists vs per-run directories: expected / completed / failed / timed-out / stalled / interrupted (progress file, no record, no live process) / running / pending.
Also lists live worker processes so orphans are visible.  python inventory.py -> records/inventory.json, records/INVENTORY.md"""
import glob, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import campaign as C
live = set()
for p in os.listdir("/proc"):
    if p.isdigit():
        try: c = open(f"/proc/{p}/cmdline").read().replace("\0", " ")
        except Exception: continue
        if c.startswith("/usr/bin/python3 -u") and "campaign.py worker" in c: live.add(c.split("worker ")[-1].strip())
inv = {}; detail = []
for f in sorted(glob.glob(os.path.join(HERE, "cases", "*.json"))):
    name = os.path.basename(f)[:-5]
    if name.endswith(".sha256"): continue
    cases = json.load(open(f)); cnt = {}
    for c in cases:
        rid = C.spec_id({k: v for k, v in c.items() if k != "commit"}); rd = os.path.join(C.RUNS, rid)
        if os.path.exists(os.path.join(rd, "record.json")): st = json.load(open(os.path.join(rd, "record.json")))["status"]
        elif rd in live: st = "running"
        elif os.path.exists(os.path.join(rd, "progress.json")) or os.path.isdir(rd): st = "interrupted"
        else: st = "pending"
        cnt[st] = cnt.get(st, 0) + 1
        if st not in ("completed", "pending"): detail.append(dict(set=name, id=rid, status=st))
    inv[name] = dict(expected=len(cases), **cnt)
all_ids = {os.path.basename(d) for d in glob.glob(os.path.join(C.RUNS, "*"))}
exp_ids = {C.spec_id({k: v for k, v in c.items() if k != "commit"}) for f in glob.glob(os.path.join(HERE, "cases", "*.json")) if not f.endswith("sha256") for c in json.load(open(f))}
extra = sorted(all_ids - exp_ids)
out = dict(sets=inv, non_completed_non_pending=detail, run_dirs_not_in_any_case_list=extra, live_workers=sorted(live), total_expected_unique=len(exp_ids), total_dirs=len(all_ids))
json.dump(out, open(os.path.join(HERE, "records", "inventory.json"), "w"), indent=1)
L = ["| set | expected | completed | failed | timeout | stalled | interrupted | running | pending |", "|---|---|---|---|---|---|---|---|---|"]
for n, c in inv.items(): L.append(f"| {n} | {c['expected']} | " + " | ".join(str(c.get(k, 0)) for k in ("completed", "failed", "timeout", "stalled", "interrupted", "running", "pending")) + " |")
open(os.path.join(HERE, "records", "INVENTORY.md"), "w").write("\n".join(L) + f"\n\nLive worker processes: {len(live)}. Run dirs not in any case list: {extra}\n"); print("\n".join(L)); print(detail, extra, len(live))
