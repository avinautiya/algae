"""Reproducible figures. python -m ice_integration.make_figures"""
import json, glob
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path("ice_integration/results"); F=R/"figures"; F.mkdir(exist_ok=True)
rec=json.load(open(R/"EVALUATION_RECORD.json"))
errs=rec["nominal"]["errs"]
bases=["constant_rate_fitted","constant_h_fitted","db04_unfitted","db03_unfitted","gnielinski_smooth_unfitted","db04_scale_fitted"]
fig,ax=plt.subplots(1,3,figsize=(12,3.6),sharey=True)
for a,s in zip(ax,["P3","P5","pooled"]):
    v=[errs[f"{b}|{s}"]["MAE_m"] for b in bases]
    a.barh(range(len(bases)),v,color=["#888"]*5+["#1f77b4"]); a.set_yticks(range(len(bases)));a.set_yticklabels(bases,fontsize=7);a.set_title(f"{s} MAE (m), nominal test");a.invert_yaxis()
plt.tight_layout();plt.savefig(F/"nominal_mae.png",dpi=130);plt.close()
w=rec["candidate_wins_over_baseline_by_station"]; bs=[b for b in w]
fig,ax=plt.subplots(figsize=(7,3.4))
for i,s in enumerate(["P3","P5","pooled"]): ax.bar([j+0.27*i for j in range(len(bs))],[w[b][s] for b in bs],0.25,label=s)
ax.axhline(81,ls=":",c="k");ax.set_xticks([j+0.27 for j in range(len(bs))]);ax.set_xticklabels(bs,rotation=20,fontsize=7);ax.set_ylabel("scenarios candidate wins (of 81)");ax.legend()
plt.tight_layout();plt.savefig(F/"scenario_wins.png",dpi=130);plt.close()
fig,ax=plt.subplots(1,2,figsize=(10,3.4))
for p in sorted(glob.glob(str(R/"demo/demo_dx*.json"))):
    if p.endswith("rows.json"):continue
    d=json.load(open(p));n=Path(p).stem.replace("demo_","");rows=d["rows"]
    ax[0].plot([r["Tw"] for r in rows],label=n);ax[1].plot([r["f_toe"] for r in rows],label=n)
ax[0].set_title("mean inflow T (°C) per 6 h");ax[1].set_title("toe flow fraction");ax[1].legend(fontsize=6)
plt.tight_layout();plt.savefig(F/"demo_dynamics.png",dpi=130);plt.close()
print("ok")
