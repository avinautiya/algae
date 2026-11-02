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
    a.barh(range(len(bases)),v,color=["#888"]*5+["#1f77b4"]); a.set_yticks(range(len(bases)));a.set_yticklabels(bases,fontsize=7);a.set_title(f"{s}: MAE (m)\n(nominal scenario; n={errs[bases[0]+'|'+s]['n']} test intervals)",fontsize=9);a.set_xlabel("mean absolute error of incision (m)");a.invert_yaxis()
fig.text(0.5,-0.02,"Single nominal scenario, one calibration. No confidence intervals: n is 4 (P3) and 2 (P5); these are retrospective, inspected data. Blue = candidate (scale-fitted DB0.4); grey = baselines (first two fitted on training data only).",ha="center",fontsize=7);plt.tight_layout();plt.savefig(F/"nominal_mae.png",bbox_inches="tight",dpi=130);plt.close()
w=rec["candidate_wins_over_baseline_by_station"]; bs=[b for b in w]
fig,ax=plt.subplots(figsize=(7,3.4))
for i,s in enumerate(["P3","P5","pooled"]): ax.bar([j+0.27*i for j in range(len(bs))],[w[b][s] for b in bs],0.25,label=s)
ax.axhline(81,ls=":",c="k");ax.set_xticks([j+0.27 for j in range(len(bs))]);ax.set_xticklabels(bs,rotation=20,fontsize=7);ax.set_xlabel("baseline the candidate (db04_scale_fitted) is compared against; MAE on test intervals",fontsize=8);ax.set_ylabel("scenarios won (of 81)");ax.set_title("Win COUNTS over engineering scenarios (clock ±2 h, T ±0.05 °C, width 1.5/2/2.5 m, 3 gap imputers)\nnot probabilities and not confidence intervals; scenarios are not independent draws; zero-height bar = 0 wins",fontsize=8);ax.legend(loc="upper right",fontsize=7)
plt.tight_layout();plt.savefig(F/"scenario_wins.png",dpi=130);plt.close()
fig,ax=plt.subplots(1,2,figsize=(10,3.4))
for p in sorted(glob.glob(str(R/"demo/demo_dx*.json"))):
    if p.endswith("rows.json"):continue
    d=json.load(open(p));n=Path(p).stem.replace("demo_","");rows=d["rows"]
    ls={"strict":"-","reconstructed":"--"}[n.split("_")[1]];ax[0].plot([r["Tw"] for r in rows],ls,label=n,lw=1.2);ax[1].plot([r["f_toe"] for r in rows],ls,label=n,lw=1.2)
ax[0].set_title("mean inflow temperature per 6 h interval");ax[0].set_ylabel("°C (supported time only)");ax[1].set_title("toe flow fraction Q_toe/Q");ax[1].set_ylabel("fraction (-)")
for a in ax: a.set_xlabel("interval k (6 h; from 2019-07-20 00:00, naive clock)")
ax[1].legend(fontsize=6);fig.text(0.5,-0.03,"Demonstration on an imposed-Q strip deposit; NOT field validation. Hydraulics were capped at 60 s per interval (not all intervals fully settled; see CONVERGENCE report). Strict/reconstructed overlap except where gaps were bridged; coarse (0.05 m) and fine (0.025 m) differ in the early contact phase.",ha="center",fontsize=7)
plt.tight_layout();plt.savefig(F/"demo_dynamics.png",bbox_inches="tight",dpi=130);plt.close()
print("ok")
