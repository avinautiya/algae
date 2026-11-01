"""Raw-trace re-scoring with the pinned external repair package (m3_comparison_repair @ 0bee4db5, imported separately; NOT modified):
  * classification.recovery.classify_run on ACTUAL raw run exports (sustained-contact labels, one-day persistence after the actual export time);
  * history_key on actual causal inputs (frozen-spec physics/boundary/geometry + deposit volume + interventions) -> shared pre-loss histories;
  * my R-based recovery-strength label (metrics.classify) kept SEPARATE;
  * refusals of incompatible exports are preserved verbatim.
`python -m ice_cliff_survival.reactivation.rescore_m3 PKG_DIR SPEC.json OUT.json DIR1 [DIR2 ...]`"""
import glob, json, os, sys
import numpy as np
def main(pkg, spec_path, out, dirs):
    sys.path.insert(0, os.path.join(pkg, "classification")); import recovery as R
    from . import metrics as M
    spec = json.load(open(spec_path)); geom = {k: v for k, v in spec["geometry"].items() if k != "dx"}; res = []
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "*.json"))):
            b = os.path.basename(f)
            if b.startswith(("queue", "launch", "jobs", "boundary", "fine_warm", "analysis")) or "_zb" in b: continue
            try: run = json.load(open(f))
            except Exception as e: res.append(dict(file=f, status="REFUSED", reason="unreadable: " + str(e))); continue
            if not isinstance(run, dict) or "records" not in run: res.append(dict(file=f, status="REFUSED", reason="not a run export")); continue
            row = dict(file=os.path.relpath(f), name=run.get("name") or b[:-5], layout=run.get("layout"), dx=run.get("dx"), removal_day=run.get("removal_day"), days=run.get("days"), incision=run.get("incision", "full"))
            try: c = R.classify_run(run)
            except (ValueError, KeyError, TypeError) as e: row.update(status="REFUSED", reason=str(e)); res.append(row); continue
            row.update(status="classified", sustained_contact_label=c["classification"], loss_day=(c["loss"] or {}).get("estimate_day") if c.get("loss") else None, loss_bracket=(c["loss"] or {}).get("bracket_day") if c.get("loss") else None, actual_removal_day=c["actual_removal_day"], removal_delay_days=c["removal_delay_days"],
                       onset=[(x["start_estimate_day"], x["onset_order_relative_to_removal"], x["bracket_confirmed_persistence_days"]) for x in c["sustained_recovery_candidates"]], truncated=c["truncated"], observed_end_day=c["observed_end_day"], preintervention_rewet=len(c["preintervention_rewet_intervals"]))
            # separate R-based strength label (needs the bare share of the same resolution/version; supplied by the caller through f_ref in the bare run of the SAME directory)
            row["_f"] = ([x["day"] for x in run["records"]], [x["f_toe"] for x in run["records"]]); row["_mean_f"] = float(np.mean(row["_f"][1])); row["_rec"] = run
            ev = row["loss_day"]; event = ev if ev is not None else float(run.get("days", row["observed_end_day"]))
            interventions = [dict(time_day=run["export"][0]["t"] / 86400, kind="clear_deposit", volume_m3=run["export"][0]["volume_m3"])] if run.get("export") else []
            row["history_key"] = R.history_key(initial_state=dict(layout=run["layout"], geometry=geom, deposit_volume_m3=run.get("deposit_volume", 0.0), incision=run.get("incision", "full")), physics=spec["physics"], boundary=spec["boundary"], forcing_history=[dict(time_day=0.0, Q=spec["physics"]["Q_in"], T=spec["physics"]["T_in"])], interventions=interventions, event_day=event)
            res.append(row)
    # R labels: bare reference per (dir, dx)
    bare = {}
    for r in res:
        if r["status"] == "classified" and r["layout"] == "bare" and r["removal_day"] is None: bare[(os.path.dirname(r["file"]), r["dx"])] = r["_mean_f"]
    for r in res:
        if r["status"] != "classified": continue
        ref = bare.get((os.path.dirname(r["file"]), r["dx"]))
        if ref and r["removal_day"] is not None: m = M.classify(np.array(r["_f"][0]), np.array(r["_f"][1]), r["actual_removal_day"], f_ref=ref); r["recovery_strength_R"] = dict(cls=m["cls"], R=m["R"])
        elif ref: r["recovery_strength_R"] = dict(cls=M.classify(np.array(r["_f"][0]), np.array(r["_f"][1]), None, f_ref=ref)["cls"], R=None)
        del r["_f"], r["_mean_f"], r["_rec"]
    groups = {}
    for r in res:
        if r["status"] == "classified" and r["loss_day"] is not None: groups.setdefault(r["history_key"], []).append(r["file"])
    out_ = dict(package_commit="0bee4db523407166cd51c85d4ed046d769864985", spec=spec_path, rows=res, shared_pre_loss_histories=[dict(key=k[:12], runs=v) for k, v in groups.items()], note="history_key ignores grid and future clearing by design; solver versions are NOT in the key (same physical input history); a hash is not proof of statistical independence")
    json.dump(out_, open(out, "w"), indent=1, default=float); return out_
if __name__ == "__main__":
    o = main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
    for r in o["rows"]: print("%-34s %-10s %-26s loss %-7s rem %-6s trunc %-5s | R-class %s" % (r.get("name", "")[:34], r["status"], r.get("sustained_contact_label", r.get("reason", "")[:26]), None if r.get("loss_day") is None else round(r["loss_day"], 3), r.get("actual_removal_day"), r.get("truncated"), (r.get("recovery_strength_R") or {}).get("cls")))
    print("shared pre-loss histories:", [(g["key"], len(g["runs"])) for g in o["shared_pre_loss_histories"]])
