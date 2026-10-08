#!/usr/bin/env python
"""
Computed vs measured Fe(III)-induced absorption of purpurogallin (tier D mechanistic check).

Tier D of the forcing model is built from the MEASURED Fe-purpurogallin spectrum (Prochazkova et al.
2025, Fig. 4; phase2/tddft_calibration.py). This script asks whether TD-DFT on the Level 3 Fe(III)
complex reproduces that increment, i.e. whether the computed photophysics supports the empirical tier:

  computed increment   d_eps(l) = eps_complex(l) - eps_ligand(l)       (same functional/basis/solvent;
                       ligand = Level 1 purpurogallin; per mole of purpurogallin)
  measured increment   dA(l) = A_PG+Fe(l) - A_PG(l)                    (same PG concentration)

Concentration-free metrics (the measured absorbances are in arbitrary path length x concentration):
  vis_ratio   mean increment over 500-750 nm / ligand peak in 280-400 nm   (computed vs measured)
  shape_r     Pearson r of the increments over 400-750 nm
  lmct_nm     wavelength of the strongest computed excitation above 450 nm with f > 0.005, with its f

    python phase1/fe_increment.py --ligand results/level1 --complex results/level3_catecholate \
        results/level3_tropolonate --outdir results/fe_increment
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
MEASURED = os.path.join(HERE, "..", "data", "empirical", "prochazkova2025_fig4_PG_PGFe_absorbance_digitized.csv")
GRID = np.arange(280.0, 800.0 + 0.1, 1.0)


def eps_of(folder, functional):
    """Molar absorptivity (L mol^-1 cm^-1) on GRID, the excited states, and the run summary."""
    level = json.load(open(os.path.join(folder, "summary.json"))).get("level") or os.path.basename(folder.rstrip("/"))
    stem = next(os.path.join(folder, f[:-len("_spectrum.csv")]) for f in sorted(os.listdir(folder))
                if f.endswith(f"_{functional}_spectrum.csv"))
    sp = pd.read_csv(stem + "_spectrum.csv").sort_values("Wavelength_nm")
    st = pd.read_csv(stem + "_states.csv")
    return np.interp(GRID, sp.Wavelength_nm, sp["Epsilon_L_mol-1_cm-1"], left=np.nan, right=np.nan), st, level


def metrics(inc, lig, band=(500, 750), peak=(280, 400), shape=(400, 750)):
    b = (GRID >= band[0]) & (GRID <= band[1])
    pk = (GRID >= peak[0]) & (GRID <= peak[1])
    return dict(vis_ratio=float(np.nanmean(inc[b]) / np.nanmax(lig[pk])))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ligand", default=os.path.join(HERE, "results", "level1"))
    p.add_argument("--complex", nargs="+", default=[os.path.join(HERE, "results", "level3_catecholate"),
                                                    os.path.join(HERE, "results", "level3_tropolonate")])
    p.add_argument("--functional", default="B3LYP")
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "fe_increment"))
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)

    m = pd.read_csv(MEASURED).sort_values("wavelength_nm")
    A_pg = np.interp(GRID, m.wavelength_nm, m.A_PG, left=np.nan)
    dA = np.interp(GRID, m.wavelength_nm, m.A_PG_Fe, left=np.nan) - A_pg
    shape = (GRID >= 400) & (GRID <= 750)
    rows = [dict(model="measured (Prochazkova et al. 2025, Fig. 4)", **metrics(dA, A_pg), shape_r=1.0,
                 lmct_nm=np.nan, lmct_f=np.nan)]
    lig, _, _ = eps_of(a.ligand, a.functional)
    curves = {"measured dA (scaled)": dA / np.nanmax(A_pg[(GRID >= 280) & (GRID <= 400)])}
    for c in a.complex:
        if not os.path.isfile(os.path.join(c, "summary.json")):
            print(f"skip {c}: no summary.json (run not finished)")
            continue
        eps_c, st, level = eps_of(c, a.functional)
        inc = eps_c - lig
        ok = shape & np.isfinite(inc) & np.isfinite(dA)
        vis = st[(st.Wavelength_nm > 450) & (st.Oscillator_Strength > 0.005)]
        top = vis.loc[vis.Oscillator_Strength.idxmax()] if len(vis) else None
        rows.append(dict(model=f"TD-{a.functional} {level}", **metrics(inc, lig),
                         shape_r=float(np.corrcoef(inc[ok], dA[ok])[0, 1]),
                         lmct_nm=np.nan if top is None else float(top.Wavelength_nm),
                         lmct_f=np.nan if top is None else float(top.Oscillator_Strength),
                         n_states_vis=int(len(vis))))
        curves[f"computed {level}"] = inc / np.nanmax(lig[(GRID >= 280) & (GRID <= 400)])
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(a.outdir, "fe_increment_metrics.csv"), index=False, float_format="%.4g")
    pd.DataFrame(dict(wavelength_nm=GRID, **curves)).to_csv(os.path.join(a.outdir, "fe_increment_curves.csv"),
                                                            index=False, float_format="%.5g")
    print(tab.round(3).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    for k, v in curves.items():
        ax.plot(GRID, v, lw=1.6 if k.startswith("measured") else 1.0, color="k" if k.startswith("measured") else None,
                label=k)
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlim(300, 800)
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("Fe-induced change / ligand UV peak")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "FigS_fe_increment.png"), dpi=600)
    fig.savefig(os.path.join(a.outdir, "FigS_fe_increment.pdf"))


if __name__ == "__main__":
    main()
