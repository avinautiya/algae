"""Station-level error diagnosis + training-only identifiability. No new fitting on any test/evaluation outcome.
Calibration is the frozen V5 one (training intervals only). Test observations enter ONLY descriptive error tables (inspected retrospective data).
python -m ice_integration.diagnosis.station_diagnosis OUT.json"""
import json, math, sys, os
from .. import point_model as PM
V5=os.path.join(os.path.dirname(__file__),"..","..","ice_swiss_field_validation","conditional_benchmark","v5_local_geometry_development")
NOM="clock+0h_T+0.00C_W2.00m"
def load(imp): return json.load(open(os.path.join(V5,imp,"CONDITIONAL_RESULTS.json")))["scenarios"]
def main(out):
    S={i:load(i) for i in("linear","left","right")}; nom=S["linear"][NOM]; cal=nom["calibration"]
    res={"test_intervals":[],"training":{},"exclusions_nominal":[]}
    # ---- per-interval descriptive table (nominal) with input-range across all completed scenarios x imputers
    for it in nom["accepted_intervals"]:
        p=PM.predict(it,cal); ts=it["temperature_support"]
        preds=[];
        for imp,sc in S.items():
            for n,s in sc.items():
                if s.get("status")!="completed": continue
                for j in s["accepted_intervals"]:
                    if j["id"]==it["id"]: preds.append(PM.predict(j,s["calibration"])["db04_scale_fitted"])
        rec=dict(id=it["id"],station=it["station"],role=it["role"],days=it["duration_s"]/86400,observed_m=it["observed_lowering_m"],obs_sigma_m=it["summed_epoch_error_scenario_m"],
          pred_m=p["db04_scale_fitted"],error_m=p["db04_scale_fitted"]-it["observed_lowering_m"],temp_observed_fraction=ts["observed_fraction"],temp_reconstructed_fraction=ts["reconstructed_fraction"],
          engineering_pred_min_m=min(preds),engineering_pred_max_m=max(preds),n_scenarios=len(preds),
          error_sign_same_in_all=all((q-it["observed_lowering_m"])*(p["db04_scale_fitted"]-it["observed_lowering_m"])>0 for q in preds),
          mean_depth_m=sum(c["local_depth_m"]*c["integral_C_s"] for c in it["hydraulic_cells"])/max(sum(c["integral_C_s"] for c in it["hydraulic_cells"]),1e-12),
          mean_Re=sum(c["Re"]*c["integral_C_s"] for c in it["hydraulic_cells"])/max(sum(c["integral_C_s"] for c in it["hydraulic_cells"]),1e-12))
        # error outside the engineering range AND outside obs sigma => not explained by clock/T/width/imputer
        lo,hi=rec["engineering_pred_min_m"],rec["engineering_pred_max_m"]; o=rec["observed_m"]; s=rec["obs_sigma_m"]
        rec["observation_inside_engineering_range_plus_sigma"]=(lo-s)<=o<=(hi+s)
        res["test_intervals"].append(rec)
    # ---- training-only identifiability (db04 scale; implied scale per training interval)
    tr=[i for i in nom["accepted_intervals"] if i["role"]=="train"]
    for i in tr:
        pass
    b=[i["db04_basis_m"] for i in tr]; y=[i["observed_lowering_m"] for i in tr]; sg=[i["summed_epoch_error_scenario_m"] for i in tr]
    implied=[yy/bb for yy,bb in zip(y,b)]
    # least-squares through origin (scale) and its sampling uncertainty from stated survey sigma (propagated, not statistical CI)
    w=lambda idx:[1.0 for _ in idx]
    def ls(idx): return sum(b[k]*y[k] for k in idx)/sum(b[k]**2 for k in idx)
    allidx=list(range(len(tr))); s_all=ls(allidx)
    loo={tr[k]["id"]:ls([j for j in allidx if j!=k]) for k in allidx}
    contrib={tr[k]["id"]:b[k]**2/sum(bb**2 for bb in b) for k in allidx}
    s_sigma=math.sqrt(sum((b[k]*sg[k])**2 for k in allidx))/sum(bb**2 for bb in b)
    by_st={}
    for st in("P3","P5"):
        idx=[k for k in allidx if tr[k]["station"]==st]; by_st[st]=dict(n_train=len(idx),ids=[tr[k]["id"] for k in idx],station_only_scale=(ls(idx) if idx else None),
            implied_scale_per_interval={tr[k]["id"]:implied[k] for k in idx},obs_sigma_m={tr[k]["id"]:sg[k] for k in idx},observed_m={tr[k]["id"]:y[k] for k in idx})
    res["training"]=dict(frozen_scale=cal["db04_scale"],recomputed_ls_scale=s_all,ls_matches_frozen=abs(s_all-cal["db04_scale"])<1e-6*max(1,s_all),
        leverage_share_of_sum_b2=contrib,leave_one_out_scale=loo,loo_range=[min(loo.values()),max(loo.values())],propagated_survey_sigma_of_scale=s_sigma,
        implied_scale_range=[min(implied),max(implied)],by_station=by_st,
        note="Frozen calibration was reproduced from training rows only. P5 contributes one 0.76-day interval whose observed lowering (0.13 m) is of the order of its 0.2 m summed survey error, so a P5-specific parameter is not identifiable from training data; no P5 fit to evaluation outcomes was attempted.")
    res["exclusions_nominal"]=[dict(id=e["id"],station=e["station"],role=e["role"],reason=e["reason"]["excluded"]) for e in nom["excluded_intervals"]]
    json.dump(res,open(out,"w"),indent=1,default=float)
    for r in res["test_intervals"]: print(r["id"],r["role"],"obs %.2f pred %.2f err %+.2f recon %.3f eng[%.2f,%.2f] inside=%s"%(r["observed_m"],r["pred_m"],r["error_m"],r["temp_reconstructed_fraction"],r["engineering_pred_min_m"],r["engineering_pred_max_m"],r["observation_inside_engineering_range_plus_sigma"]))
    t=res["training"]; print({k:t[k] for k in("frozen_scale","recomputed_ls_scale","ls_matches_frozen","loo_range","propagated_survey_sigma_of_scale","implied_scale_range")}); print(t["leverage_share_of_sum_b2"]); print(t["by_station"])
if __name__=="__main__": main(sys.argv[1])
