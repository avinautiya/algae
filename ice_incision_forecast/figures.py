"""Reproduce six method figures from a saved synthetic forecast; no CFD results."""
import argparse
import csv
import hashlib
import itertools
import json
import os
import platform
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"ice-incision-mpl-cache"))
os.environ.setdefault("XDG_CACHE_HOME",str(Path(tempfile.gettempdir())/"ice-incision-xdg-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .core import ConductancePrior,ThermalMap,region_bounds
from .forecast import forecast_case,profile_powers,validate


def csv_file(path,rows):
    rows=list(rows)
    with path.open("x",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def integrated_profile(case,profile):
    heat=np.zeros(len(profile))
    for f in case["forcing"]:
        model=ThermalMap(f["inlet_Q_m3_s"],f["inlet_excess_K"],**case["properties"])
        heat+=np.asarray(profile_powers(model,profile,f["conductance_multiplier"]))*f["duration_s"]
    areas=np.asarray([c["width_m"]*c["length_m"] for c in case["cells"]])
    return heat/(case["properties"]["rho_ice"]*case["properties"]["latent_heat"]*areas)


def regional_interval(case,prior,start,stop):
    w=region_bounds(prior,start,stop); values=[]
    area=sum(c["width_m"]*c["length_m"] for c in case["cells"][start:stop])
    for label in ("minimum","maximum"):
        energy=0.0
        for f in case["forcing"]:
            model=ThermalMap(f["inlet_Q_m3_s"],f["inlet_excess_K"],**case["properties"])
            multiplier=f["conductance_multiplier"]
            energy+=model.region_power(multiplier*w[label]["prefix_W_K"],multiplier*w[label]["region_W_K"])*f["duration_s"]
        values.append(1000*energy/(case["properties"]["rho_ice"]*case["properties"]["latent_heat"]*area))
    return values


def independent_rk4(model,g,steps):
    p=model.inlet_power_W; q=model.inlet_Q_m3_s; step=g/steps
    def rhs(state):
        p,q=state; v=p/(model.rho_water*model.cp*q)
        return np.asarray([-v,v/(model.rho_water*model.latent_heat)])
    state=np.asarray([p,q])
    for _ in range(steps):
        k1=rhs(state); k2=rhs(state+step*k1/2); k3=rhs(state+step*k2/2); k4=rhs(state+step*k3)
        state+=step*(k1+2*k2+2*k3+k4)/6
    return state


def make_figures(results,output):
    case=json.loads((results/"case.json").read_text())
    recorded=json.loads((results/"manifest.json").read_text())
    digest=hashlib.sha256(json.dumps(case,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
    if digest!=recorded["case_sha256"]:
        raise ValueError("Case differs from the saved forecast manifest")
    if output.exists():
        raise FileExistsError("Choose a new figure directory")
    prior=validate(case); result=forecast_case(case)
    if len(case["regions"])!=1:
        raise ValueError("This six-figure demonstration expects one declared region")
    if prior.exact_prefix_W_K or prior.total_lower_W_K!=prior.total_upper_W_K:
        raise ValueError("The six-figure demonstration requires an exact-total no-prefix case")
    if not 0<case["regions"][0]["start"]<case["regions"][0]["stop"]<len(case["cells"]):
        raise ValueError("The six-figure demonstration requires an interior target region")
    output.mkdir(parents=True); data=output/"data"; data.mkdir()
    plt.rcParams.update({"font.size":9,"axes.labelsize":10,"axes.titlesize":11,"legend.fontsize":8,
        "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"savefig.dpi":300})
    figures=[]
    def save(fig,name,caption):
        fig.text(0.5,0.008,"Synthetic conditional-model demonstration",ha="center",fontsize=8,color="0.35")
        fig.tight_layout(rect=(0,0.025,1,1))
        for suffix in ("png","pdf"):
            fig.savefig(output/(name+"."+suffix),bbox_inches="tight")
        plt.close(fig); figures.append(dict(name=name,caption=caption,formats=["png","pdf"]))
    lengths=np.asarray([c["length_m"] for c in case["cells"]]); edges=np.r_[0,np.cumsum(lengths)]
    x=(edges[:-1]+edges[1:])/2; region=case["regions"][0]; a,b=region["start"],region["stop"]
    w=result["regional_witnesses"][region["name"]]
    f=case["forcing"][0]; model=ThermalMap(f["inlet_Q_m3_s"],f["inlet_excess_K"],**case["properties"])
    fig,axes=plt.subplots(3,1,figsize=(6.7,7.6),sharex=True)
    rows=[]
    for label,color in (("minimum","#3266a8"),("maximum","#c44c35")):
        profile=np.asarray(w[label]["profile_W_K"])
        axes[0].step(edges,np.r_[profile,profile[-1]],where="post",label=label+" regional melt",color=color)
        power=[model.state(g*f["conductance_multiplier"])["power_W"]/1000 for g in np.r_[0,np.cumsum(profile)]]
        axes[1].plot(edges,power,color=color,label=label+" witness")
        lowering=integrated_profile(case,profile)*1000
        axes[2].step(edges,np.r_[lowering,lowering[-1]],where="post",color=color)
        for i in range(len(profile)):
            rows.append(dict(witness=label,cell=i,x_mid_m=x[i],conductance_W_K=profile[i],lowering_1h_mm=lowering[i]))
    for ax in axes:
        ax.axvspan(edges[a],edges[b],color="0.7",alpha=.2)
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Cell conductance (W/K)"); axes[0].legend(loc="upper right")
    axes[1].set_ylabel("Sensible power (kW)"); axes[1].set_title("Same total conductance gives identical outlet states")
    axes[2].set_ylabel("Bed lowering (mm)"); axes[2].set_xlabel("Downstream distance (m)")
    axes[0].set_title("Attaining profiles for the declared control region")
    csv_file(data/"fig01_witness_profiles.csv",rows)
    save(fig,"fig01_equal_outlet_different_incision","Two admissible profiles attain opposite regional extrema while sharing outlet states. Shading marks the predeclared region. Power is for the first interval; lowering integrates the full horizon.")

    last=[r for r in result["bed"] if r["time_s"]==result["summary"]["horizon_s"]]
    lo=np.asarray([r["lowering_min_m"] for r in last])*1000
    hi=np.asarray([r["lowering_max_m"] for r in last])*1000
    nominal=np.asarray([r["lowering_nominal_m"] for r in last])*1000
    fig,ax=plt.subplots(figsize=(6.7,4.2))
    ax.fill_between(x,lo,hi,color="#3266a8",alpha=.2,label="Sharp marginal range")
    rng=np.random.default_rng(20261031); profile=np.asarray([c["conductance_nominal_W_K"] for c in case["cells"]])
    sampled=[]
    # Illustration for the no-prefix exact-total demonstration only.
    for step in range(2500):
        i,j=rng.choice(len(profile),2,replace=False)
        low=max(prior.lower_W_K[i]-profile[i],profile[j]-prior.upper_W_K[j])
        high=min(prior.upper_W_K[i]-profile[i],profile[j]-prior.lower_W_K[j])
        delta=rng.uniform(low,high); profile[i]+=delta; profile[j]-=delta
        if step>=1000 and step%10==0:
            y=integrated_profile(case,profile)*1000
            sampled.append(y.tolist()); ax.plot(x,y,color="#3266a8",alpha=.045,lw=.6)
    ax.plot(x,nominal,color="#1b2733",lw=1.8,label="Nominal profile")
    ax.axvspan(edges[a],edges[b],color="0.7",alpha=.18)
    ax.set(xlabel="Downstream distance (m)",ylabel="Bed lowering (mm)",title="Local forecast uncertainty survives exact total conductance")
    ax.legend(); ax.grid(alpha=.2)
    csv_file(data/"fig02_marginal_envelope.csv",[dict(cell=i,x_m=x[i],lower_mm=lo[i],nominal_mm=nominal[i],upper_mm=hi[i]) for i in range(len(x))])
    (data/"fig02_sampled_profiles.json").write_text(json.dumps(dict(seed=20261031,lowering_mm=sampled,
        note="Feasible-profile illustration, not probability samples or confidence intervals"))+"\n")
    save(fig,"fig02_sharp_spatial_envelope","Pointwise bounds are marginal and cannot all be attained together. Sampled feasible profiles illustrate the class without defining a probability distribution.")

    fig,axes=plt.subplots(2,1,figsize=(6.7,5.8),sharex=True)
    times=np.r_[0,np.cumsum([r["duration_s"] for r in case["forcing"]])]/3600
    axes[0].step(times,np.r_[[r["inlet_Q_m3_s"] for r in case["forcing"]],case["forcing"][-1]["inlet_Q_m3_s"]],where="post",color="#3266a8")
    twin=axes[0].twinx(); twin.step(times,np.r_[[r["inlet_excess_K"] for r in case["forcing"]],case["forcing"][-1]["inlet_excess_K"]],where="post",color="#c44c35")
    axes[0].set_ylabel("Inlet discharge (m³/s)",color="#3266a8"); twin.set_ylabel("Inlet excess (K)",color="#c44c35")
    rrows=[r for r in result["regions"] if r["region"]==region["name"]]
    rlo=np.r_[0,[r["lowering_min_m"]*1000 for r in rrows]]
    rhi=np.r_[0,[r["lowering_max_m"]*1000 for r in rrows]]
    rnom=np.r_[0,[r["lowering_nominal_m"]*1000 for r in rrows]]
    rmid=np.r_[0,[r["lowering_minimax_m"]*1000 for r in rrows]]
    axes[1].fill_between(times,rlo,rhi,alpha=.22,color="#3266a8",label="Sharp regional range")
    axes[1].plot(times,rnom,color="0.5",ls="--",label="Nominal conductance profile")
    axes[1].plot(times,rmid,color="#1b2733",marker="o",label="Scalar minimax forecast")
    if "critical_lowering_m" in region:
        axes[1].axhline(region["critical_lowering_m"]*1000,color="#c44c35",ls=":",label="Declared incision threshold")
    axes[1].set(xlabel="Forecast time (h)",ylabel="Regional mean lowering (mm)")
    axes[0].set_title("A forward forecast driven by prescribed inlet scenarios")
    axes[1].legend(); axes[1].grid(alpha=.2)
    csv_file(data/"fig03_regional_forecast.csv",[dict(time_h=t,lower_mm=l,nominal_mm=m,minimax_mm=c,upper_mm=u) for t,l,m,c,u in zip(times,rlo,rnom,rmid,rhi)])
    save(fig,"fig03_regional_forecast","A single fixed witness attains each time-integrated extreme despite changing inlet scenarios and the common conductance multiplier. The midpoint minimizes scalar worst-case absolute error in this class. The declared incision threshold is not a flood threshold; no future outlet measurement is an input.")

    sensor=[]; base=regional_interval(case,prior,a,b)
    true_g=[c["conductance_nominal_W_K"] for c in case["cells"]]
    for k in range(1,len(true_g)):
        p=ConductancePrior(prior.lower_W_K,prior.upper_W_K,prior.total_lower_W_K,prior.total_upper_W_K,((k,sum(true_g[:k])),))
        lower,upper=regional_interval(case,p,a,b)
        sensor.append(dict(position_m=edges[k],width_mm=upper-lower))
    bracket=ConductancePrior(prior.lower_W_K,prior.upper_W_K,prior.total_lower_W_K,prior.total_upper_W_K,
        ((a,sum(true_g[:a])),(b,sum(true_g[:b]))))
    bracket_bounds=regional_interval(case,bracket,a,b)
    fig,ax=plt.subplots(figsize=(6.7,4.0))
    ax.axhline(base[1]-base[0],color="0.4",ls="--",label="Outlet information only")
    ax.plot([r["position_m"] for r in sensor],[r["width_mm"] for r in sensor],marker="o",ms=3,color="#3266a8",label="One ideal interior prefix observation")
    ax.axhline(max(0,bracket_bounds[1]-bracket_bounds[0]),color="#c44c35",ls=":",lw=2,label="Both region boundaries observed exactly")
    ax.axvspan(edges[a],edges[b],color="0.7",alpha=.18)
    ax.set(xlabel="Ideal interior observation position (m)",ylabel="Regional forecast range width (mm)",title="Which thermal information resolves local incision?")
    ax.legend(loc="lower right"); ax.grid(alpha=.2)
    csv_file(data/"fig04_ideal_prefix_information.csv",sensor)
    save(fig,"fig04_thermal_information_limit","Idealized exact prefix information partitions the uncertainty class. Real point sensors and noisy measurements require observation operators and uncertainty treatment; this is not a field sensor-performance claim.")

    verification=ThermalMap(.02,.3); g=100000; exact=verification.state(g)["power_W"]
    ns=np.asarray([8,16,32,64]); errors=np.asarray([abs(independent_rk4(verification,g,int(n))[0]-exact)/verification.inlet_power_W for n in ns])
    fig,axes=plt.subplots(1,3,figsize=(9.0,3.2))
    axes[0].loglog(ns,errors,"o-",color="#3266a8"); axes[0].loglog(ns,errors[0]*(ns/ns[0])**(-4),"--",color="0.4",label="Fourth-order guide")
    axes[0].set(xlabel="RK4 steps",ylabel="Relative power error",title="Independent ODE check"); axes[0].legend()
    axes[0].set_xticks(ns); axes[0].set_xticklabels([str(int(n)) for n in ns]); axes[0].minorticks_off()
    p=ConductancePrior((0,)*4,(30000,)*4,50000,50000)
    candidates=[[10000*x for x in v] for v in itertools.product(range(4),repeat=4) if sum(v)==5]
    pairs=[]
    for start in range(4):
        for stop in range(start+1,5):
            rb=region_bounds(p,start,stop)
            vals=[verification.region_power(sum(v[:start]),sum(v[start:stop])) for v in candidates]
            for label,val in (("minimum",min(vals)),("maximum",max(vals))):
                wv=rb[label]; bound=verification.region_power(wv["prefix_W_K"],wv["region_W_K"])
                pairs.append(dict(start=start,stop=stop,kind=label,exhaustive_W=val,theorem_W=bound))
    axes[1].scatter([r["exhaustive_W"]/1000 for r in pairs],[r["theorem_W"]/1000 for r in pairs],s=16,color="#3266a8")
    lim=max(r["theorem_W"] for r in pairs)/1000; axes[1].plot([0,lim],[0,lim],"--",color="0.4")
    axes[1].set(xlabel="Enumerated extreme (kW)",ylabel="Sharp formula (kW)",title="Attainment check")
    stes=np.geomspace(1e-5,.02,35); corrections=[]
    for ste in stes:
        m=ThermalMap(.02,ste*334000/4217); gg=1.2*m.capacity_W_K
        corrections.append(100*(m.state(gg)["power_W"]/m.inlet_power_W-np.exp(-1.2))/np.exp(-1.2))
    axes[2].semilogx(stes,corrections,color="#c44c35")
    axes[2].set(xlabel="Sensible/latent Stefan number",ylabel="Power correction (%)",title="Meltwater mass addition")
    for ax in axes: ax.grid(alpha=.2)
    csv_file(data/"fig05_rk4_convergence.csv",[dict(steps=int(n),relative_error=e) for n,e in zip(ns,errors)])
    csv_file(data/"fig05_exhaustive_extrema.csv",pairs)
    csv_file(data/"fig05_mass_addition.csv",[dict(Ste=ste,power_correction_percent=cor) for ste,cor in zip(stes,corrections)])
    save(fig,"fig05_independent_checks","Independent RK4 integration, exhaustive small-case extrema and the physical effect of meltwater addition check different parts of the implementation. Numerical agreement supports the implementation; the proof establishes the conditional theorem.")

    local=np.linspace(.05,.8,16); total=np.linspace(0,.5,12); records=[]; images=[]
    broad_start=max(1,len(true_g)//10); broad_stop=len(true_g)-broad_start
    targets=[("Declared control",a,b),("Broader reach average",broad_start,broad_stop)]
    for label,start,stop in targets:
        image=np.zeros((len(local),len(total)))
        for i,spread in enumerate(local):
            for j,band in enumerate(total):
                p=ConductancePrior(tuple(g*(1-spread) for g in true_g),tuple(g*(1+spread) for g in true_g),
                    sum(true_g)*(1-band),sum(true_g)*(1+band))
                lower,upper=regional_interval(case,p,start,stop); image[i,j]=upper-lower
                records.append(dict(target=label,start=start,stop=stop,local_conductance_halfwidth=spread,
                    total_conductance_halfwidth=band,regional_range_mm=upper-lower))
        images.append(image)
    fig=plt.figure(figsize=(9,4.1)); grid=fig.add_gridspec(1,3,width_ratios=[1,1,.06])
    axes=[fig.add_subplot(grid[0,j]) for j in range(2)]
    limit=max(float(np.max(im)) for im in images)
    for ax,image,target in zip(axes,images,targets):
        im=ax.imshow(image,origin="lower",aspect="auto",extent=(total[0],total[-1],local[0],local[-1]),cmap="viridis",vmin=0,vmax=limit)
        ax.set(xlabel="Total conductance half-width",ylabel="Local conductance half-width",title=target[0])
    cb=fig.colorbar(im,cax=fig.add_subplot(grid[0,2])); cb.set_label("Forecast range width (mm)")
    csv_file(data/"fig06_uncertainty_map.csv",records)
    save(fig,"fig06_local_vs_total_information","Information sensitivity for the declared control and a secondary broader-reach observable. Bulk information can constrain a reach average while leaving local allocation uncertain. The second region is an explanatory scale comparison, not a replacement primary endpoint.")
    manifest=dict(case_sha256=digest,figures=figures,scope="All figures use the supplied synthetic reduced-model case; no field or CFD results",
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),matplotlib_version=matplotlib.__version__,
        numpy_version=np.__version__,python_version=platform.python_version())
    (output/"figure_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    (output/"CAPTIONS.md").write_text("# Reproducible method figures\n\nAll panels are synthetic conditional-model demonstrations.\n\n"+"\n\n".join("**"+r["name"]+".** "+r["caption"] for r in figures)+"\n")
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path,required=True); parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    try:
        m=make_figures(args.results,args.output)
        print(json.dumps(dict(figures=len(m["figures"]),output=str(args.output)),indent=2))
    except (ValueError,KeyError,OSError) as exc:
        parser.exit(2,f"Figure generation failed: {exc}\n")


if __name__=="__main__":
    main()
