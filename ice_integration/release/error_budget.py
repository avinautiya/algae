"""Signed error-budget numbers + calibration-influence figure. Reads frozen results only; no refit to evaluation intervals.
python -m ice_integration.release.error_budget OUTDIR"""
import json, math, os, sys, re
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..",".."))
REC=json.load(open(os.path.join(ROOT,"ice_integration","results","EVALUATION_RECORD.json")))
DIAG=json.load(open(os.path.join(ROOT,"ice_integration","results","STATION_DIAGNOSIS.json")))
C="db04_scale_fitted"
def factor_effects():
    rows=[]
    for key,v in REC["table"].items():
        if v.get("status")!="completed": continue
        imp,name=key.split("|"); m=re.match(r"clock([+-]\d+)h_T([+-][\d.]+)C_W([\d.]+)m",name)
        rows.append(dict(imp=imp,clock=int(m[1]),T=float(m[2]),W=float(m[3]),pooled_bias=v["errs"][f"{C}|pooled"]["bias_m"],pooled_mae=v["errs"][f"{C}|pooled"]["MAE_m"],P3_bias=(v["errs"].get(f"{C}|P3") or {}).get("bias_m"),P5_bias=(v["errs"].get(f"{C}|P5") or {}).get("bias_m")))
    out={}
    for f in("imp","clock","T","W"):
        g=defaultdict(list)
        for r in rows: g[r[f]].append(r)
        out[f]={str(k):dict(n=len(v),mean_pooled_bias_m=sum(x["pooled_bias"] for x in v)/len(v),min_pooled_bias_m=min(x["pooled_bias"] for x in v),max_pooled_bias_m=max(x["pooled_bias"] for x in v),mean_pooled_MAE_m=sum(x["pooled_mae"] for x in v)/len(v)) for k,v in sorted(g.items())}
    out["all"]=dict(n=len(rows),pooled_bias_range_m=[min(r["pooled_bias"] for r in rows),max(r["pooled_bias"] for r in rows)],fraction_scenarios_with_positive_pooled_bias=sum(r["pooled_bias"]>0 for r in rows)/len(rows))
    return out
def depth_sensitivity(Q=1.2,w=2.0,T=1.0):
    mu=1.79e-3; rho=1000.0; k=0.56  # water near 0 C
    def h_ht(d):
        v=Q/(w*d); Dh=4*w*d/(w+2*d); Re=rho*v*Dh/mu; Nu=0.023*Re**0.8*(13.0**0.4); return k*Nu/Dh
    base=0.62; return {f"depth {d:+.2f} m vs proxy {base} m":h_ht(base+d)/h_ht(base)-1 for d in(-0.1,0.1,0.2,0.3,0.4)}
def influence(outdir):
    tr=[i for i in DIAG["test_intervals"] if i["role"]=="train"]; T=DIAG["training"]
    ids=[i["id"] for i in tr]; w=[T["leverage_share_of_sum_b2"][i] for i in ids]; res=[i["error_m"] for i in tr]; loo=[T["leave_one_out_scale"][i] for i in ids]
    fig,ax=plt.subplots(1,3,figsize=(12,3.6))
    ax[0].bar(ids,res,color=["#1f77b4" if s=="P3" else "#d62728" for s in [i["station"] for i in tr]]);ax[0].axhline(0,c="k",lw=.6);ax[0].set_ylabel("training residual, pred − obs (m)");ax[0].set_title("Training residuals (frozen scale 4.40)\nblue P3, red P5; survey σ ±0.2 m not drawn");ax[0].tick_params(axis="x",rotation=45,labelsize=8)
    ax[1].bar(ids,[x*100 for x in w],color="#555");ax[1].set_ylabel("share of Σ(basis²) (%)");ax[1].set_title("Calibration weight per interval\n(scale ∝ Σ b·y / Σ b²)");ax[1].tick_params(axis="x",rotation=45,labelsize=8)
    ax[2].bar(ids,loo,color="#888");ax[2].axhline(T["frozen_scale"],c="k",ls="--",lw=1,label="frozen scale 4.40 (all six)");ax[2].set_ylabel("DB0.4 scale estimate (–)");ax[2].set_title("Leave-one-interval-out scale\n(descriptive range, not a confidence interval)");ax[2].legend(fontsize=7);ax[2].tick_params(axis="x",rotation=45,labelsize=8)
    fig.text(.5,-.04,"Training rows only. No evaluation interval enters any panel. Six training intervals; effective sample ≈ 1.3 (P3:3 carries 85 %). P5 has one 0.76-day interval below its survey error.",ha="center",fontsize=7)
    plt.tight_layout(); plt.savefig(os.path.join(outdir,"calibration_influence.png"),dpi=130,bbox_inches="tight"); plt.close()
    return dict(ids=ids,residual_m=res,weight=w,loo_scale=loo)
def main(outdir):
    os.makedirs(outdir,exist_ok=True)
    out=dict(factor_effects=factor_effects(),depth_proxy_sensitivity_fractional_change_in_heat_transfer=depth_sensitivity(),calibration_influence=influence(outdir),
      test_signed_errors_nominal_m={i["id"]:i["error_m"] for i in DIAG["test_intervals"] if i["role"]=="test"})
    json.dump(out,open(os.path.join(outdir,"ERROR_BUDGET_NUMBERS.json"),"w"),indent=1); print(json.dumps({k:out[k] for k in("factor_effects","depth_proxy_sensitivity_fractional_change_in_heat_transfer")},indent=1)[:3500])
if __name__=="__main__": main(sys.argv[1])
