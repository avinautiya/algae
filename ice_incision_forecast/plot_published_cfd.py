"""Replot completed, pinned CFD metrics; never executes or modifies OpenFOAM."""
import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"ice-incision-mpl-cache"))
os.environ.setdefault("XDG_CACHE_HOME",str(Path(tempfile.gettempdir())/"ice-incision-xdg-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SOURCE="https://github.com/avinautiya/algae/blob/9a52524/ice_network_v2/cfd/graetz_final.json"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new directory")
    source=Path(__file__).resolve().parent/"records/published_graetz_9a52524.json"
    raw=source.read_bytes(); data=json.loads(raw); rows=data["space"]
    output=args.output; output.mkdir(parents=True)
    plt.rcParams.update({"font.size":9,"axes.labelsize":10,"axes.titlesize":11,"legend.fontsize":8,
                         "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    fig,ax=plt.subplots(figsize=(6.7,4.4))
    columns=[("err_u_over_Ub","Velocity / bulk speed","o"),
             ("err_T_shape_rms","Temperature shape / bulk excess","s"),
             ("err_wall_flux_rel","Wall heat flux / reference","^")]
    h=[r["dy_over_H"] for r in rows]; orders={}
    for key,label,marker in columns:
        values=[r[key] for r in rows]
        ax.loglog(h,values,marker=marker,label=label)
        orders[key]=[math.log(a/b)/math.log(ha/hb) for a,b,ha,hb in zip(values,values[1:],h,h[1:])]
    guide=[rows[0]["err_u_over_Ub"]*(v/h[0])**2 for v in h]
    ax.loglog(h,guide,ls="--",color="0.4",alpha=.6,label="Second-order guide")
    ax.set_xticks(h); ax.set_xticklabels(["1/10","1/20","1/40"]); ax.minorticks_off()
    ax.set(xlabel="Cross-channel cell size / channel height",ylabel="Reported normalized error",
           title="Committed OpenFOAM laminar verification")
    ax.grid(alpha=.2); ax.legend(loc="upper left")
    fig.text(.5,.01,"Repository metrics at 9a52524 — independently replotted, not rerun",ha="center",fontsize=8,color="0.35")
    fig.tight_layout(rect=(0,.03,1,1))
    for suffix in ("png","pdf"):
        fig.savefig(output/("fig07_published_cfd_verification."+suffix),dpi=300,bbox_inches="tight")
    plt.close(fig)
    base=data["time"][0]; sensitivity={}
    for key,label,_ in columns:
        absolute=max(abs(r[key]-base[key]) for r in data["time"])
        sensitivity[key]=dict(max_absolute_change_in_error_metric=absolute,
            change_over_base_spatial_error=absolute/base[key] if base[key] else None)
    report=dict(source_url=SOURCE,source_sha256=hashlib.sha256(raw).hexdigest(),observed_orders=orders,
        steady_sensitivity=sensitivity,
        note="Stored temperature and wall-flux error metrics do not all differ by less than 1e-13. Differences remain tiny relative to spatial error; this is a reporting correction, not a demonstrated physical failure.",
        scope="Replot of prior completed verification; no new CFD run or active-work modification")
    (output/"source_and_sensitivity.json").write_text(json.dumps(report,indent=2)+"\n")
    with (output/"plotted_metrics.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    (output/"CAPTION.md").write_text("# Committed CFD verification\n\nSpatial refinement of reported velocity, temperature-shape and wall-flux errors from the completed laminar benchmark at commit `9a52524`. Each curve uses its own stated normalization. This is an independent replot of stored metrics, not an independent solver rerun. Source URL, hash, observed orders and steady-state sensitivity are in `source_and_sensitivity.json`.\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
