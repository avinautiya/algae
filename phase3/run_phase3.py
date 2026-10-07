#!/usr/bin/env python
"""
Phase 3 driver: Monte Carlo uncertainty propagation (Latin Hypercube) and global
sensitivity analysis (Saltelli / Sobol') through the Phase 2 chain
(TD-DFT spectrum -> packaging -> cell optics -> BioSNICAR -> forcing).

Examples
--------
python run_phase3.py --phase1-l2 ../phase1/results/level2             # full analysis
python run_phase3.py --demo --n-mc 400 --n-sobol 256                  # quick pipeline check
python run_phase3.py --phase1-l2 ... --reuse                           # re-plot / re-analyse cached runs

Cost (17 ms per model run on one CPU core; --workers parallelises):
  Monte Carlo        n_mc                          (1000)
  Sobol' (2nd order) n_sobol * (2D + 2)            (1024 * 24 = 24,576 for D = 11)
  Sobol' by scale    n_sobol * (G + 2), G = 3      (5,120)
  Tornado            2D + 1                        (23)
  -> about 10 min on 1 core, ~5 min on a 2-vCPU Colab runtime.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from param_space import ParameterSpace, default_parameters  # noqa: E402
import stats_tools as st  # noqa: E402

ALL_OUTPUTS = ["rf_A", "rf_B", "rf_C", "rf_D", "eff_A", "eff_B", "eff_C", "eff_D",
               "bba_A", "bba_B", "bba_C", "bba_D", "d_CB", "d_DC", "d_BA", "d_CA", "d_DA"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--functional", default="B3LYP")
    p.add_argument("--demo", action="store_true", help="BioSNICAR ppg.csv instead of Phase 1 (labelled DEMO)")
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--outdir", default=os.path.join(HERE, "results"))
    p.add_argument("--n-mc", type=int, default=1000, help="Latin Hypercube sample size")
    p.add_argument("--n-sobol", type=int, default=1024, help="Saltelli base sample (power of 2)")
    p.add_argument("--no-second-order", action="store_true")
    p.add_argument("--no-groups", action="store_true", help="skip the scale-grouped Sobol' design")
    p.add_argument("--no-tier-d-params", action="store_true", help="drop the Fe-complexed fraction parameter")
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--seed", type=int, default=2024)
    p.add_argument("--sza", type=float, default=45.0, help="deg; ~solar noon at S6 in July")
    p.add_argument("--sw-down", type=float, default=None)
    p.add_argument("--sobol-outputs", nargs="+",
                   default=["rf_A", "rf_B", "rf_C", "rf_D", "eff_C", "eff_D", "d_CB", "d_DC"])
    p.add_argument("--fig3b-outputs", nargs=2, default=["rf_C", "eff_C"])
    p.add_argument("--group-outputs", nargs="+", default=["rf_C", "rf_D", "eff_C", "eff_D"])
    p.add_argument("--tornado-output", default="rf_D")
    p.add_argument("--surface-output", default="eff_C")
    p.add_argument("--s2-output", default="eff_C")
    p.add_argument("--reuse", action="store_true", help="reuse cached model runs in outdir/tables")
    p.add_argument("--usetex", action="store_true")
    return p.parse_args(argv)


def run_or_load(path, design, model_kwargs, workers, reuse, label):
    from forward_model import evaluate_many
    if reuse and os.path.isfile(path):
        df = pd.read_csv(path)
        if len(df) == len(design):
            print(f"{label}: reusing {path}")
            return df
        print(f"{label}: cached file has {len(df)} rows, need {len(design)} -> recomputing")
    print(f"{label}: {len(design)} model runs")
    df = evaluate_many(model_kwargs, design, workers=workers)
    df.to_csv(path, index=False, float_format="%.8g")
    return df


def main(argv=None):
    a = parse_args(argv)
    t0 = time.time()
    tab = os.path.join(a.outdir, "tables")
    figd = os.path.join(a.outdir, "figures")
    os.makedirs(tab, exist_ok=True)
    os.makedirs(figd, exist_ok=True)

    # empirical calibration of the Phase 1 spectrum -> molecular PDFs (phase2/tddft_calibration.py)
    import biosnicar_bridge as bb
    import cell_optics as co
    import tddft_calibration as TC
    lig = co.demo_spectrum(bb.locate_biosnicar(a.biosnicar)) if a.demo else co.load_phase1(a.phase1_l2, "level2",
                                                                                          a.functional)
    cal = TC.calibrate(lig)
    cal_summary = cal.summary()
    with open(os.path.join(tab, "tddft_calibration.json"), "w") as fh:
        json.dump(cal_summary, fh, indent=1)
    space = ParameterSpace(default_parameters(include_tier_d=not a.no_tier_d_params, calibration=cal_summary))
    space.table().to_csv(os.path.join(tab, "parameters.csv"), index=False)
    print(f"{space.D} uncertain parameters: {', '.join(space.names)}")
    model_kwargs = dict(phase1_l2=None if a.demo else a.phase1_l2, functional=a.functional, demo=a.demo,
                        biosnicar=a.biosnicar, sza=a.sza, sw_down=a.sw_down,
                        qtable_cache=os.path.join(tab, "qstar_table.npz"), calibration_point=cal.point())
    second = not a.no_second_order

    # ---------------------------------------------------------- 1) Monte Carlo (LHS)
    X = space.lhs(a.n_mc, seed=a.seed)
    mc = run_or_load(os.path.join(tab, "mc_lhs_runs.csv"), space.as_dicts(X), model_kwargs,
                     a.workers, a.reuse, "Monte Carlo (LHS)")
    summary = st.describe_frame(mc, [c for c in ALL_OUTPUTS if c in mc])
    summary.to_csv(os.path.join(tab, "mc_summary.csv"), float_format="%.6g")
    print("\nMonte Carlo summary (W m^-2 for rf/d/eff; albedo for bba):")
    print(summary.loc[["rf_A", "rf_B", "rf_C", "rf_D", "d_CB", "d_DC"],
                      ["mean", "median", "p2_5", "p97_5", "q25", "q75", "mean_ci_lo", "mean_ci_hi"]]
          .round(2).to_string())

    # ---------------------------------------------------------- 2) Sobol' (parameters)
    Xs, prob = space.saltelli(a.n_sobol, second_order=second, seed=a.seed)
    sb = run_or_load(os.path.join(tab, "saltelli_runs.csv"), space.as_dicts(Xs), model_kwargs,
                     a.workers, a.reuse, "Sobol' (parameter level)")
    sobol_tabs, s2_tabs, dominant = {}, {}, []
    for out in a.sobol_outputs:
        tabo, s2 = st.sobol(prob, sb[out].to_numpy(), second_order=second, seed=a.seed)
        tabo.to_csv(os.path.join(tab, f"sobol_{out}.csv"), float_format="%.5f")
        if s2 is not None:
            s2.to_csv(os.path.join(tab, f"sobol_S2_{out}.csv"), float_format="%.5f")
        sobol_tabs[out], s2_tabs[out] = tabo, s2
        top = tabo.ST.idxmax()
        dominant.append(dict(output=out, top_parameter=top, top_ST=tabo.ST.max(), top_S1=tabo.S1[top],
                             sum_S1=tabo.S1.sum(), second_parameter=tabo.ST.nlargest(2).index[-1],
                             second_ST=tabo.ST.nlargest(2).iloc[-1]))
    dom = pd.DataFrame(dominant)
    dom.to_csv(os.path.join(tab, "sobol_dominant_parameters.csv"), index=False, float_format="%.4f")
    print("\nDominant parameters (total-effect index S_T):")
    print(dom.round(3).to_string(index=False))
    conv = st.sobol_convergence(prob, sb[a.fig3b_outputs[0]].to_numpy(), a.n_sobol, second)
    conv.to_csv(os.path.join(tab, f"sobol_convergence_{a.fig3b_outputs[0]}.csv"), index=False)

    # ---------------------------------------------------------- 3) Sobol' (scales)
    group_tabs = {}
    if not a.no_groups:
        Xg, probg = space.saltelli(a.n_sobol, second_order=False, groups=True, seed=a.seed + 1)
        sg = run_or_load(os.path.join(tab, "saltelli_group_runs.csv"), space.as_dicts(Xg), model_kwargs,
                         a.workers, a.reuse, "Sobol' (scale groups)")
        rows = []
        for out in a.group_outputs:
            tg, _ = st.sobol(probg, sg[out].to_numpy(), second_order=False, seed=a.seed)
            group_tabs[out] = tg
            for g, r in tg.iterrows():
                rows.append(dict(output=out, group=g, **r.to_dict()))
        gdf = pd.DataFrame(rows)
        gdf.to_csv(os.path.join(tab, "sobol_groups.csv"), index=False, float_format="%.5f")
        print("\nScale-level total effects S_T (molecular / cellular / environmental):")
        print(gdf.pivot(index="output", columns="group", values="ST").round(3).to_string())

    # ---------------------------------------------------------- 4) Tornado (OAT)
    from forward_model import evaluate_many
    td = evaluate_many(model_kwargs, st.tornado_design(space), workers=1, progress=False)
    torn = st.tornado_table(td, a.tornado_output)
    torn.to_csv(os.path.join(tab, f"tornado_{a.tornado_output}.csv"), index=False, float_format="%.5f")

    # ---------------------------------------------------------- 5) Figures
    import figures3 as F3
    F3.set_style(a.usetex)
    F3.save(F3.fig3a(mc, demo=a.demo), figd, "Fig3A_forcing_distributions")
    if group_tabs:
        F3.save(F3.fig3b({k: sobol_tabs[k] for k in a.fig3b_outputs}, group_tabs, space, demo=a.demo),
                figd, "Fig3B_sobol_indices")
    F3.save(F3.fig3c(torn, mc, s2_tabs.get(a.s2_output), space, a.tornado_output, a.surface_output,
                     a.s2_output, demo=a.demo), figd, "Fig3C_interactions")
    F3.save(F3.fig_s3_convergence(conv, space, a.fig3b_outputs[0]), figd, "FigS3_sobol_convergence")

    cfg = vars(a).copy()
    cfg.update(parameters=space.table().drop(columns="rationale").to_dict(orient="records"),
               n_model_runs=int(len(mc) + len(sb) + (len(Xg) if not a.no_groups else 0) + len(td)),
               runtime_s=round(time.time() - t0, 1))
    with open(os.path.join(a.outdir, "run_config.json"), "w") as fh:
        json.dump(cfg, fh, indent=2, default=str)
    print(f"\nDone in {time.time() - t0:.0f} s -> {a.outdir}")


if __name__ == "__main__":
    main()
