"""Paper-ready plots from the predeclared synthetic study; no CFD launched."""
import argparse
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR","/tmp/ice-thermal-observability-mpl")
os.environ.setdefault("XDG_CACHE_HOME","/tmp/ice-thermal-observability-cache")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from .channel import Channel,CrossStream
from .forecast import dump
from .recovery import sharp_interval

BLUE="#28668c"; RED="#be503b"; TEAL="#208879"; DARK="#263540"


def make(study,output):
    study=Path(study);out=Path(output)
    if out.exists():raise FileExistsError("Choose a new figure directory")
    def read(name):return json.loads((study/name).read_text())
    protocol=read("protocol.json");summary=read("summary.json")
    out.mkdir(parents=True);(out/"data").mkdir()
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                         "axes.spines.right":False,"axes.titleweight":"semibold",
                         "axes.labelcolor":DARK,"text.color":DARK,"axes.grid":True,
                         "grid.alpha":.16,"figure.facecolor":"white","savefig.facecolor":"white"})
    saved=[]
    def save(fig,name,data):
        if len(fig.axes)>1:
            for i,ax in enumerate(fig.axes):
                ax.text(0,1.04,"("+chr(97+i)+")",transform=ax.transAxes,fontsize=10,fontweight="bold",va="bottom")
        fig.text(.5,.008,"SYNTHETIC FIXED-FLOW MODEL",ha="center",fontsize=8,color="#666666")
        fig.tight_layout(rect=[0,.035,1,.94])
        fig.savefig(out/(name+".png"),dpi=180)
        fig.savefig(out/(name+".pdf"))
        dump(out/"data"/(name+".json"),data)
        plt.close(fig);saved.append(name)
    channel=Channel();model=CrossStream(channel,protocol["primary_cells"])
    stations=protocol["stations_m"];region=protocol["target_m"]
    cap=protocol["temperature_cap_K"]
    forecast=read("forecast/forecast.json")["regions"][0]
    # Figure 1: analytical, smooth, positive continuum counterexample.
    witness=read("continuum_witness.json");pc=Channel(velocity="plug")
    coeff_lo=np.array(witness["minus_coefficients_K"]);coeff_hi=np.array(witness["plus_coefficients_K"])
    count=len(coeff_lo);wave=(np.arange(count)+.5)*np.pi/pc.depth_m
    lam=pc.diffusivity*wave**2/pc.speed_m_s
    y=np.linspace(0,pc.depth_m,501);x=np.linspace(0,.16,501)
    theta_lo=np.sin(np.outer(y,wave))@coeff_lo;theta_hi=np.sin(np.outer(y,wave))@coeff_hi
    decay=np.exp(-np.outer(x,lam));bulk_delta=decay@((coeff_hi-coeff_lo)/(pc.depth_m*wave))
    flux_delta=pc.conductivity*decay@(wave*(coeff_hi-coeff_lo))
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    fig.suptitle("Smooth positive profiles with identical bulk calibration",fontsize=15)
    axes[0].plot(theta_lo,y*1000,color=BLUE,label="State A")
    axes[0].plot(theta_hi,y*1000,color=RED,ls="--",label="State B")
    axes[0].set(xlabel="Upstream temperature excess (K)",ylabel="Height above ice (mm)")
    axes[0].legend(frameon=False)
    axes[1].plot(x,bulk_delta*1000,color=TEAL)
    axes[1].axhline(0,color="#777777",lw=.8)
    axes[1].scatter(stations,[0,0],color=DARK,zorder=5,label="Identical observations")
    axes[1].set(xlabel="Distance (m)",ylabel="Bulk temperature B − A (mK)")
    axes[1].legend(frameon=False,fontsize=8)
    axes[2].plot(x,flux_delta,color=RED)
    axes[2].axhspan(0,0,color="none")
    axes[2].axvspan(*region,color=BLUE,alpha=.12,label="Forecast reach")
    axes[2].axhline(0,color="#777777",lw=.8)
    axes[2].set(xlabel="Distance (m)",ylabel="Bed heat flux B − A (W m$^{-2}$)")
    axes[2].legend(frameon=False,fontsize=8)
    save(fig,"fig01_continuum_counterexample",dict(y_m=y.tolist(),state_A_K=theta_lo.tolist(),state_B_K=theta_hi.tolist(),
         x_m=x.tolist(),bulk_difference_K=bulk_delta.tolist(),bed_flux_difference_W_m2=flux_delta.tolist(),witness=witness))
    # Figure 2: attaining discrete shear-flow profiles, propagated to current section.
    low=np.array(forecast["lower_witness"]["initial_profile_K"])
    high=np.array(forecast["upper_witness"]["initial_profile_K"])
    probe=max(stations);a,b=region
    tlo=model.propagate(low,probe);thi=model.propagate(high,probe)
    xx=np.linspace(.012,.15,241)
    flo=np.array([model.wall_flux_W_m2(low,z) for z in xx])
    fhi=np.array([model.wall_flux_W_m2(high,z) for z in xx])
    factor=model.lowering_factor(a,b,protocol["horizon_s"])*1e6
    interval=forecast["power_interval"]
    fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    fig.suptitle("Sharp regional erosion forecasts under no-slip laminar shear",fontsize=15)
    axes[0].plot(tlo,model.y*1000,color=BLUE,label="Minimum-erosion witness")
    axes[0].plot(thi,model.y*1000,color=RED,ls="--",label="Maximum-erosion witness")
    axes[0].set(xlabel="Current temperature excess (K)",ylabel="Height above ice (mm)")
    axes[0].legend(frameon=False,fontsize=8,loc="lower right")
    axes[1].plot(xx,flo,color=BLUE);axes[1].plot(xx,fhi,color=RED,ls="--")
    axes[1].axvspan(a,b,color=TEAL,alpha=.12)
    axes[1].set(xlabel="Distance (m)",ylabel="Bed heat flux (W m$^{-2}$)")
    lo_um=interval["lower"]*factor;hi_um=interval["upper"]*factor;mid=interval["minimax"]*factor
    axes[2].plot([0,1],[lo_um,lo_um],color=BLUE,lw=2)
    axes[2].plot([0,1],[hi_um,hi_um],color=RED,lw=2)
    axes[2].fill_between([0,1],lo_um,hi_um,color=TEAL,alpha=.15)
    axes[2].scatter([.5],[mid],color=DARK,s=65,zorder=4,label="Minimax forecast")
    axes[2].set(xlim=(-.1,1.1),xticks=[],ylabel="60 s regional lowering (µm)")
    axes[2].text(.5,.1,f"{lo_um:.3f}–{hi_um:.3f} µm",transform=axes[2].transAxes,ha="center")
    axes[2].text(.5,.64,"Minimax forecast",transform=axes[2].transAxes,ha="center",fontsize=9)
    save(fig,"fig02_sharp_shear_forecast",dict(y_m=model.y.tolist(),current_lower_witness_K=tlo.tolist(),
         current_upper_witness_K=thi.tolist(),x_m=xx.tolist(),lower_witness_flux_W_m2=flo.tolist(),
         upper_witness_flux_W_m2=fhi.tolist(),lowering_range_um=[lo_um,hi_um],minimax_um=mid))
    # Figure 3: deterministic measurement error and target-specific profile information.
    noise=read("measurement_noise.json");c=model.heat_kernel_W_K(a,b)
    ec=np.array([model.bulk_kernel(s) for s in stations]);reference=np.full(model.cells,.1)
    current_kernel=model.heat_kernel_W_K(a-probe,b-probe)
    # c_orig = S(probe)^T c_current; this moment acts on the current profile.
    goal_row=c/current_kernel.sum();goal=[]
    for sigma in protocol["measurement_noise_K"]:
        ee=np.vstack([ec,goal_row]);dd=ee@reference;eps=np.r_[np.zeros(2),sigma]
        rr=sharp_interval(c,ee,dd,eps,np.zeros(model.cells),np.full(model.cells,cap))
        goal.append(dict(noise_bound_K=sigma,radius_W=rr["worst_case_absolute_error"]))
    fig,axes=plt.subplots(1,2,figsize=(10,3.8))
    fig.suptitle("Information requirements for a regional ice-erosion forecast",fontsize=15)
    positions=np.arange(len(protocol["measurement_noise_K"]))
    for j,count in enumerate((1,2)):
        rows=[r for r in noise if r["stations"]==count]
        radii=np.array([(r["upper_W"]-r["lower_W"])/2 for r in rows])*factor
        axes[0].plot(positions,radii,"o-",color=(BLUE,RED)[j],label=f"{count} bulk station"+("s" if count>1 else ""))
    axes[0].set(xticks=positions,xticklabels=["0","0.01","0.1","1"],xlabel="Absolute temperature-error bound (mK)",ylabel="Conditional lowering error radius (µm)")
    axes[0].legend(frameon=False,fontsize=9)
    base_radius=interval["worst_case_absolute_error"]*factor
    axes[1].axhline(base_radius,color=RED,ls="--",label="Two exact bulk observations")
    axes[1].plot(positions,[g["radius_W"]*factor for g in goal],"o-",color=TEAL,label="Plus target-weighted profile moment")
    axes[1].set(xticks=positions,xticklabels=["0","0.01","0.1","1"],xlabel="Profile-moment error bound (mK)",ylabel="Conditional lowering error radius (µm)")
    axes[1].legend(frameon=False,fontsize=8,loc="center left")
    save(fig,"fig03_measurement_information",dict(noise=noise,goal_information=goal,
         meaning="Illustrative ideal weighted-profile observation; not one point thermometer or an achieved sensor precision"))
    # Figure 4: every predeclared baseline, no selection of a favorable score.
    cases=read("smooth_profile_cases.json");methods=["uniform_entrance","two_mode_graetz","minimax"]
    labels=["Uniform entrance","Two-mode Graetz","Conditional minimax"]
    truth=np.array([v["true_power_W"] for v in cases]);fig,axes=plt.subplots(1,3,figsize=(12,3.8))
    fig.suptitle("Comparison on 120 seeded smooth admissible inlet profiles",fontsize=15)
    predictions=[]
    for ax,name,label,color in zip(axes,methods,labels,(BLUE,RED,TEAL)):
        pred=np.array([v["predictions_W"][name] for v in cases]);predictions.append(pred)
        bounds=[min(truth.min(),pred.min()),max(truth.max(),pred.max())]
        ax.scatter(truth,pred,s=13,color=color,alpha=.65)
        ax.plot(bounds,bounds,color="#888888",lw=.8)
        mae=summary["smooth_profile_metrics"][name]["mean_abs_error_W"]
        ax.set(title=label,xlabel="Resolved-model regional heat (W)",ylabel="Predicted regional heat (W)")
        ax.text(.04,.95,f"MAE = {mae*1000:.3f} mW",transform=ax.transAxes,va="top",fontsize=9)
    save(fig,"fig04_all_baselines",dict(truth_W=truth.tolist(),predictions_W={name:v.tolist() for name,v in zip(methods,predictions)},
         summary=summary["smooth_profile_metrics"],meaning="Conditional minimax optimizes worst-case class error, not mean error under this synthetic profile distribution"))
    # Figure 5: independent analytical reference, mesh convergence, LP dual checks.
    refin=read("refinement.json");checks=read("independent_vertex_checks.json")
    fig,axes=plt.subplots(1,3,figsize=(12,3.8));fig.suptitle("Independent verification of thermal propagation and sharp bounds",fontsize=15)
    plug=[r for r in refin if r["velocity"]=="plug"]
    nn=np.array([r["cells"] for r in plug]);err=np.array([r["relative_uniform_power_error"] for r in plug])
    axes[0].loglog(nn,err,"o-",color=BLUE,label="FV versus analytical sine series")
    axes[0].loglog(nn,err[0]*(nn[0]/nn)**2,"--",color="#999999",label="Second-order guide")
    axes[0].set(xlabel="Transverse cells",ylabel="Relative regional-power error",xticks=nn,xticklabels=nn)
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].legend(frameon=False,fontsize=8)
    fixed_data=read("fixed_data_refinement.json") if (study/"fixed_data_refinement.json").exists() else None
    shear=[r for r in (fixed_data["rows"] if fixed_data else refin) if r["velocity"]=="shear_free"]
    axes[1].plot([r["cells"] for r in shear],[r["lower_W"] for r in shear],"o-",color=BLUE,label="Sharp lower bound")
    axes[1].plot([r["cells"] for r in shear],[r["upper_W"] for r in shear],"o-",color=RED,label="Sharp upper bound")
    axes[1].set(xlabel="Transverse cells",ylabel="Shear-flow regional heat (W)",title="Fixed observation data" if fixed_data else "Grid-generated data")
    axes[1].legend(frameon=False,fontsize=8)
    ids=[r["case"] for r in checks]
    axes[2].semilogy(ids,[max(r["max_error"],1e-18) for r in checks],".",color=TEAL,label="Independent vertex error")
    axes[2].semilogy(ids,[max(abs(r["max_dual_gap"]),1e-18) for r in checks],"+",color=RED,ms=4,label="Evaluated primal–dual gap")
    axes[2].set(xlabel="Seeded continuous test case",ylabel="Absolute scalar-functional error")
    axes[2].legend(frameon=False,fontsize=8)
    axes[2].text(.02,.04,"Zeros plotted at 10$^{-18}$",transform=axes[2].transAxes,fontsize=8,bbox=dict(facecolor="white",alpha=.85,edgecolor="none"))
    save(fig,"fig05_independent_verification",dict(refinement=refin,fixed_data_refinement=fixed_data,independent_vertex_checks=checks))
    # Figure 6: uncertainty changes with thermal development and observation density.
    age=read("thermal_age.json");fig,ax=plt.subplots(figsize=(7.5,4.2))
    fig.suptitle("Finite upstream information and thermal-development scale",fontsize=15)
    for count,color in ((1,BLUE),(2,RED),(4,TEAL)):
        rows=[r for r in age if r["stations"]==count]
        ax.semilogx([r["beta"] for r in rows],[100*r["half_width_over_midpoint"] for r in rows],"o-",color=color,label=f"{count} exact bulk station"+("s" if count>1 else ""))
    ax.set(xlabel="Spatial scale β = x$_1$ / (UH$^2$/α)",ylabel="Conditional half-width / midpoint (%)")
    ax.legend(frameon=False)
    if any(r.get("numerically_unresolved_width",False) for r in age):
        ax.text(.42,.65,"Some four-station ranges are below numerical resolution",transform=ax.transAxes,fontsize=8)
    save(fig,"fig06_thermal_age",dict(rows=age,definition="Probe=2β UH²/α; target=[3β,6β] UH²/α; station sets nested and upstream"))
    # Figure 7: nonnegative source responses and conservative forecast extension.
    source=read("surface_heating_forecast/forecast.json")["regions"][0]
    starts=np.linspace(.04,.14,81);source_length=.01
    fractions=[model.surface_heat_coefficient(s,s+source_length,a,b)/(channel.width_m*source_length) for s in starts]
    fig,axes=plt.subplots(1,2,figsize=(10,3.8));fig.suptitle("Surface heating reaches the ice through resolved transverse transport",fontsize=15)
    axes[0].plot(starts,100*np.array(fractions),color=RED)
    axes[0].axvspan(a,b,color=BLUE,alpha=.12,label="Target reach")
    axes[0].set(xlabel="Start of a 1 cm surface-heating patch (m)",ylabel="Patch heat delivered to target bed (%)")
    axes[0].legend(frameon=False,fontsize=8)
    for i,(data,label,color) in enumerate(((forecast,"No extra surface heat",BLUE),(source,"0–100 W m$^{-2}$ surface heat",RED))):
        vals=data["lowering_m"];center=vals["minimax"]*1e6;radius=vals["worst_case_absolute_error"]*1e6
        axes[1].errorbar(i,center,yerr=radius,fmt="o",capsize=6,color=color)
    axes[1].set(xticks=[0,1],xticklabels=["No extra heating","Bounded surface heating"],ylabel="60 s regional lowering (µm)")
    axes[1].tick_params(axis="x",labelsize=8)
    save(fig,"fig07_surface_heat",dict(source_start_m=starts.tolist(),source_length_m=source_length,
         target_heat_fraction=fractions,base_forecast=forecast["lowering_m"],heated_forecast=source["lowering_m"]))
    captions="""# Figure captions

All seven figures are synthetic computations in a fixed prescribed-flow heat equation. They are not field observations or new OpenFOAM runs. Source data are in `data/`; every plot has PNG and vector PDF exports.

1. **Continuum counterexample.** Analytical smooth positive plug-flow inlet profiles share the exact bulk readings at 0.02 and 0.04 m. Their downstream bed heat differs. The explicit finite sine expansion is independent of the finite-volume/LP implementation. This demonstrates existence, not a large uncertainty for all measurements.
2. **Sharp shear-flow forecast.** Two attaining finite-volume inlet witnesses are propagated to the current section, where they are smooth and satisfy the same bulk calibration. They attain different regional heat/lowering over 0.06–0.12 m. The bounds are conditional on the inlet box and fixed-flow model. The shaded scalar range is deterministic and jointly attained at its endpoints.
3. **Measurement information.** Left: one/two upstream bulk observations with bounded errors. Right: an additional ideal target-weighted current-profile moment, derived from the adjoint kernel. This illustrative sensing calculation is not an achieved point-sensor precision and requires transverse profile information.
4. **All baseline results.** All 120 seeded smooth inlet profiles and all predeclared comparators are retained. Both conventional comparators use the same bulk data as the interval method. Average accuracy and worst-case conditional guarantees are different criteria. A modal reconstruction may violate the admitted inlet bounds; the count is in the study summary.
5. **Independent checks.** Analytical plug-flow series, shear-grid refinement with the same observation vector on every mesh (supplemental fixed-data audit), and independent active-face enumeration/evaluated LP duals. The shear data use a 768-cell reference with the 384-to-768 change reported separately as an indicator, not a certified error bound. Zero scalar errors are displayed at a stated plotting floor. Grid convergence is not a certified continuum or environmental error enclosure.
6. **Thermal age.** Conditional relative range under nested one/two/four station sets. Physical lengths scale with UH²/alpha. Ranges flagged below numerical resolution do not establish exact identification. This is a synthetic scale study with fixed thermal and hydraulic assumptions, not a measured environmental regime boundary.
7. **Surface heating.** The fraction of prescribed surface heat reaching the target bed depends on where it enters. Source power in the reporting reach is included in the conservative heat ledger. The interval extension admits independent nonnegative surface heating of 0–100 W/m² from 0.04 to 0.12 m. It is not a radiation/evaporation/weather model.
"""
    (out/"CAPTIONS.md").write_text(captions)
    manifest=dict(study_summary_sha256=hashlib.sha256((study/"summary.json").read_bytes()).hexdigest(),
                  protocol_sha256=hashlib.sha256((study/"protocol.json").read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),matplotlib=matplotlib.__version__,numpy=np.__version__,
                  figures=saved,outputs_sha256={str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest()
                                               for f in sorted(out.rglob("*")) if f.is_file()})
    dump(out/"manifest.json",manifest)
    print(json.dumps(dict(output=str(out),figures=len(saved)),indent=2))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--study",required=True);p.add_argument("--output",required=True)
    a=p.parse_args();make(a.study,a.output)
