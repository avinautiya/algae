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
  strongest_vis_state_nm / _f   wavelength and f of the strongest computed excitation above 450 nm
              (f > 0.005). Its character (LMCT, d-d, ligand pi-pi*) is NOT assigned here - that needs a
              transition-density / NTO analysis - so the old name 'lmct_nm' was withdrawn.

Spectra are rebuilt from the stick lists (states.csv) with ONE broadening on ONE shared grid (not
interpolated from each run's own pre-broadened spectrum, whose grid may not cover 280 nm). The runs must
share functional, basis, solvent, eps, TDA/RPA and width (checked; mismatch raises), and must carry a
valid completion marker (provisional geometry allowed only with --allow-provisional; xTB geometries are
labelled exploratory).

Caveats:
  * concentration: the measured increment is per mole of purpurogallin at the measured Fe:PG ratio; if
    complexation was incomplete the measured dA is a fraction of the full-complex increment, so vis_ratio
    compares a lower bound with the computed full-complex value. shape_r is concentration-free.
  * this is a mechanistic investigation; tier D of the forcing model uses the MEASURED increment only, and
    nothing here feeds tier D.

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


SETTINGS = ("basis", "solvent", "eps", "tda", "cartesian_d")


def eps_of(folder, functional, fwhm_ev=0.3, allow_provisional=False, reference=None):
    """Molar absorptivity (L mol^-1 cm^-1) on GRID rebuilt from the sticks, the states, the level and the
    run settings. Validates the completion marker and, if `reference` settings are given, that this run
    used the same electronic-structure settings."""
    import completion
    import spectra
    ok, status, why = completion.validate_run(folder, allow_provisional=allow_provisional)
    if not ok:
        raise ValueError(f"{folder}: not a valid run ({status}: {'; '.join(why)})")
    summ = json.load(open(os.path.join(folder, "summary.json")))
    if functional not in summ.get("tddft", {}):
        raise ValueError(f"{folder}: no {functional} result (have {sorted(summ.get('tddft', {}))})")
    settings = {k: summ.get(k) for k in SETTINGS}
    if reference is not None:
        diff = {k: (settings[k], reference[k]) for k in SETTINGS if settings[k] != reference[k]}
        if diff:
            raise ValueError(f"{folder}: settings differ from the ligand run: {diff}")
    level = summ.get("level") or os.path.basename(folder.rstrip("/"))
    st = pd.read_csv(os.path.join(folder, f"{level}_{functional}_states.csv"))
    eps = spectra.gaussian_broaden(spectra.nm_to_ev(GRID), st.Energy_eV.to_numpy(float),
                                   st.Oscillator_Strength.to_numpy(float), fwhm_ev)
    cov = spectra.root_count_sensitivity(st.Energy_eV.to_numpy(float), st.Oscillator_Strength.to_numpy(float),
                                         1.0, fwhm_ev, window_nm=(GRID[0], GRID[-1]), norm_nm=(400.0, 750.0))
    return eps, st, level, dict(settings=settings, status=status, root_coverage=cov,
                                geometry=summ.get("geometry_status", {}).get("status"))


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
    p.add_argument("--fwhm", type=float, default=0.3, help="one Gaussian FWHM (eV) for ligand and complexes")
    p.add_argument("--allow-provisional", action="store_true",
                   help="accept runs whose geometry is provisional/exploratory (reported per row)")
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "fe_increment"))
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)

    m = pd.read_csv(MEASURED).sort_values("wavelength_nm")
    A_pg = np.interp(GRID, m.wavelength_nm, m.A_PG, left=np.nan)
    dA = np.interp(GRID, m.wavelength_nm, m.A_PG_Fe, left=np.nan) - A_pg
    shape = (GRID >= 400) & (GRID <= 750)
    rows = [dict(model="measured (Prochazkova et al. 2025, Fig. 4)", **metrics(dA, A_pg), shape_r=1.0,
                 strongest_vis_state_nm=np.nan, strongest_vis_state_f=np.nan)]
    lig, _, _, lig_meta = eps_of(a.ligand, a.functional, a.fwhm, a.allow_provisional)
    curves = {"measured dA (scaled)": dA / np.nanmax(A_pg[(GRID >= 280) & (GRID <= 400)])}
    for c in a.complex:
        try:
            eps_c, st, level, meta = eps_of(c, a.functional, a.fwhm, a.allow_provisional, lig_meta["settings"])
        except (OSError, ValueError, StopIteration) as e:
            print(f"skip {c}: {e}")
            continue
        inc = eps_c - lig
        ok = shape & np.isfinite(inc) & np.isfinite(dA)
        vis = st[(st.Wavelength_nm > 450) & (st.Oscillator_Strength > 0.005)]
        top = vis.loc[vis.Oscillator_Strength.idxmax()] if len(vis) else None
        rows.append(dict(model=f"TD-{a.functional} {level}", **metrics(inc, lig),
                         shape_r=float(np.corrcoef(inc[ok], dA[ok])[0, 1]),
                         strongest_vis_state_nm=np.nan if top is None else float(top.Wavelength_nm),
                         strongest_vis_state_f=np.nan if top is None else float(top.Oscillator_Strength),
                         n_states_vis=int(len(vis)), run_status=meta["status"], geometry=meta["geometry"],
                         root_cov_max_rel_change=meta["root_coverage"]["max_rel_change_in_window"],
                         exploratory=meta["geometry"] == "exploratory_xtb"))
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
