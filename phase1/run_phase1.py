#!/usr/bin/env python
"""
Phase 1 driver: geometry optimisation -> TD-DFT -> broadened MAC spectrum.

Examples
--------
# Level 1, B3LYP/6-31G(d)/IEF-PCM(water) geometry, B3LYP and CAM-B3LYP spectra:
python run_phase1.py --level level1 --tddft-functionals B3LYP CAM-B3LYP

# Level 2 on a Colab GPU runtime:
python run_phase1.py --level level2 --gpu --tddft-functionals B3LYP CAM-B3LYP

# Resume an interrupted optimisation from its last geometry:
python run_phase1.py --level level2 --gpu --start-xyz results/level2/opt_B3LYP_last.xyz

# Skip the optimisation (geometry already converged) and only redo spectra:
python run_phase1.py --level level1 --start-xyz results/level1/opt_B3LYP_final.xyz --skip-opt

Outputs in --outdir (default results/<level>):
  opt_<F>_traj.xyz / _last.xyz / _final.xyz     optimisation trajectory and result
  <level>_<F>_states.csv     [State, Wavelength_nm, Energy_eV, Oscillator_Strength, MAC_estimated, Transitions]
  <level>_<F>_spectrum.csv   [Wavelength_nm, Energy_eV, Oscillator_Strength, Epsilon_L_mol-1_cm-1, MAC_estimated]
  <level>_<F>_spectrum.npz   plot-ready NumPy arrays
  <level>_<F>_spectrum.png   MAC(lambda) with stick spectrum
  summary.json               energies, settings, timings
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import molecules
import qc
import spectra


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--level", choices=list(molecules.MOLECULES), default="level1")
    p.add_argument("--basis", default="6-31g*", help="6-31G(d) by default")
    p.add_argument("--spherical", action="store_true",
                   help="use spherical (5d) instead of Gaussian-style Cartesian (6d) d functions")
    p.add_argument("--opt-functional", default="B3LYP", help="functional for the ground-state geometry")
    p.add_argument("--tddft-functionals", nargs="+", default=["B3LYP"],
                   help="one or more functionals for vertical excitations at the optimised geometry "
                        "(e.g. B3LYP CAM-B3LYP)")
    p.add_argument("--solvent", choices=["pcm", "smd", "none"], default="pcm")
    p.add_argument("--eps", type=float, default=qc.WATER_EPS_STATIC, help="static dielectric for IEF-PCM")
    p.add_argument("--nstates", type=int, default=30)
    p.add_argument("--tda", action="store_true", help="Tamm-Dancoff approximation (faster, f less reliable)")
    p.add_argument("--fwhm", type=float, default=0.3, help="Gaussian FWHM in eV")
    p.add_argument("--lam-min", type=float, default=300.0)
    p.add_argument("--lam-max", type=float, default=800.0)
    p.add_argument("--decadic", action="store_true", help="report decadic instead of Napierian MAC")
    p.add_argument("--start-xyz", default=None, help="start (or resume) from this .xyz file (last frame)")
    p.add_argument("--rebuild", action="store_true", help="regenerate the start geometry with RDKit")
    p.add_argument("--no-xtb", action="store_true", help="skip the GFN2-xTB/ALPB pre-optimisation")
    p.add_argument("--skip-opt", action="store_true", help="no DFT optimisation (single point + TD-DFT)")
    p.add_argument("--xtb-geometry", action="store_true",
                   help="GFN2-xTB/ALPB geometry, then DFT single point + TD-DFT (no DFT optimisation)")
    p.add_argument("--maxsteps", type=int, default=100)
    p.add_argument("--freq", action="store_true", help="harmonic frequencies at the optimised geometry")
    p.add_argument("--gpu", action="store_true", help="run DFT/TD-DFT on GPU via gpu4pyscf")
    p.add_argument("--no-df", action="store_true", help="disable density fitting")
    p.add_argument("--grid-level", type=int, default=3,
                   help="DFT integration grid (3 = production; 2 is ~1.5x faster on CPU, <0.1 mEh error)")
    p.add_argument("--lebedev", type=int, default=29,
                   help="PCM Lebedev order per atomic sphere (29 -> 302 pts; 17 -> 110 pts, faster)")
    p.add_argument("--outdir", default=None)
    p.add_argument("--max-memory", type=int, default=10000,
                   help="PySCF memory limit in MB (loosely honoured: on a 16 GB machine use about 7000)")
    p.add_argument("--verbose", type=int, default=3, help="PySCF verbosity (4 prints every SCF cycle with timings)")
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    spec = molecules.MOLECULES[a.level]
    outdir = a.outdir or os.path.join("results", a.level)
    os.makedirs(outdir, exist_ok=True)
    solvent = None if a.solvent == "none" else a.solvent
    grid_opts = dict(grids_level=a.grid_level, lebedev_order=a.lebedev)
    nthreads = qc.set_threads()
    t_start = time.time()
    summary = dict(level=a.level, name=spec["name"], formula=spec["formula"],
                   molar_mass_g_mol=spec["molar_mass"], charge=spec["charge"],
                   multiplicity=spec["multiplicity"], basis=a.basis, cartesian_d=not a.spherical,
                   solvent=a.solvent, eps=a.eps if a.solvent == "pcm" else None,
                   nstates=a.nstates, tda=a.tda, fwhm_ev=a.fwhm, napierian=not a.decadic,
                   gpu=a.gpu, cpu_threads=nthreads, grid_level=a.grid_level,
                   pcm_lebedev_order=a.lebedev, tddft={})

    # ---------------- geometry ------------------------------------------------
    if a.start_xyz:
        xyz = qc.xyz_file_to_atom_block(a.start_xyz)
        print(f"Starting geometry: {a.start_xyz}")
    else:
        xyz = molecules.get_xyz(a.level, rebuild=a.rebuild)
        print("Starting geometry: built-in MMFF94s conformer" + (" (rebuilt)" if a.rebuild else ""))
    mol = molecules.get_mole(a.level, basis=a.basis, xyz=xyz, cart=not a.spherical,
                             max_memory=a.max_memory, verbose=a.verbose)
    print(f"{spec['name']}: {mol.natm} atoms, {mol.nelectron} electrons, {mol.nao} basis functions, "
          f"charge {mol.charge}, multiplicity {mol.spin + 1}", flush=True)

    if (not a.skip_opt or a.xtb_geometry) and not a.no_xtb and not a.start_xyz:
        try:
            mol = qc.xtb_preoptimize(mol, solvent="water" if solvent else None)
            with open(os.path.join(outdir, "xtb_geometry.xyz"), "w") as fh:
                fh.write(qc._mole_to_xyz(mol, f"{qc.xtb_preoptimize.method}/ALPB(water) geometry"))
            summary["xtb_method"] = qc.xtb_preoptimize.method
        except ImportError:
            print("tblite not installed -> skipping xTB pre-optimisation (pip install tblite)")
    if a.xtb_geometry:
        a.skip_opt = True
        summary["geometry"] = f"{getattr(qc.xtb_preoptimize, 'method', 'GFN2-xTB')}/ALPB(water)"

    # ---------------- ground-state optimisation -------------------------------
    t0 = time.time()
    if a.skip_opt:
        mf = qc.make_mf(mol, a.opt_functional, solvent, a.eps, density_fit=not a.no_df, use_gpu=a.gpu,
                        chkfile=os.path.join(outdir, f"scf_{a.opt_functional}.chk"), **grid_opts)
        e_gs = mf.kernel()
        if not mf.converged:
            raise RuntimeError("ground-state SCF did not converge")
    else:
        mol, mf, e_gs, opt_conv = qc.optimize_geometry(mol, a.opt_functional, solvent, a.eps,
                                             maxsteps=a.maxsteps, workdir=outdir,
                                             tag=f"opt_{a.opt_functional}", use_gpu=a.gpu,
                                             density_fit=not a.no_df, **grid_opts)
        summary["optimisation_converged"] = bool(opt_conv)
    summary["opt_functional"] = a.opt_functional
    summary["ground_state_energy_Eh"] = float(e_gs)
    summary["time_ground_state_s"] = round(time.time() - t0, 1)

    if a.freq:
        if a.gpu or not a.no_df:
            # analytic Hessian is implemented for conventional integrals on CPU
            mf_h = qc.make_mf(mol, a.opt_functional, solvent, a.eps, density_fit=False, **grid_opts)
            mf_h.kernel()
        else:
            mf_h = mf
        freqs, n_imag = qc.harmonic_check(mf_h, outdir)
        summary["n_imaginary_frequencies"] = n_imag

    # ---------------- TD-DFT + spectra for each functional --------------------
    for func in a.tddft_functionals:
        t0 = time.time()
        if func.upper() == a.opt_functional.upper():
            mf_td = mf
        else:
            mf_td = qc.make_mf(mol, func, solvent, a.eps, density_fit=not a.no_df, use_gpu=a.gpu,
                               chkfile=os.path.join(outdir, f"scf_{func}.chk"), **grid_opts)
            mf_td.kernel()
            if not mf_td.converged:
                raise RuntimeError(f"{func} SCF did not converge")
        res = qc.run_tddft(mf_td, nstates=a.nstates, tda=a.tda,
                           checkpoint=os.path.join(outdir, f"tda_{func}_x0.npz") if a.tda else None)

        lines_df, spec_df, arrays = spectra.build_spectrum(
            res["energies_ev"], res["osc_strengths"], spec["molar_mass"],
            a.lam_min, a.lam_max, step_nm=1.0, fwhm_ev=a.fwhm, napierian=not a.decadic)
        lines_df["Transitions"] = [
            "; ".join(f"{h}->{p} ({w:.2f})" for h, p, w in qc.dominant_transitions(res["td"], i))
            for i in range(len(lines_df))]

        stem = os.path.join(outdir, f"{a.level}_{func}")
        lines_df.to_csv(stem + "_states.csv", index=False, float_format="%.6g")
        spec_df.to_csv(stem + "_spectrum.csv", index=False, float_format="%.6g")
        np.savez(stem + "_spectrum.npz", **arrays)
        try:
            spectra.plot_spectrum(arrays, title=f"{spec['name'][:45]} - TD-{func}/{a.basis}, {a.solvent.upper()}",
                                  path=stem + "_spectrum.png")
        except ImportError:
            print("matplotlib not installed -> skipping PNG (CSV/NPZ written)", flush=True)

        vis = spec_df[(spec_df.Wavelength_nm >= 400) & (spec_df.Wavelength_nm <= 700)]
        i_max = int(spec_df.MAC_estimated.values.argmax())
        summary["tddft"][func] = dict(
            time_s=round(time.time() - t0, 1),
            highest_state_eV=float(res["energies_ev"].max()),
            covers_lam_min=bool(res["energies_ev"].max() + 2 * a.fwhm >= spectra.HC_EV_NM / a.lam_min),
            lambda_max_in_window_nm=float(spec_df.Wavelength_nm.iloc[i_max]),
            mac_max_m2_kg=float(spec_df.MAC_estimated.iloc[i_max]),
            mac_mean_400_700_m2_kg=float(vis.MAC_estimated.mean()),
            brightest_state=lines_df.loc[lines_df.Oscillator_Strength.idxmax()].drop("Transitions").to_dict(),
        )
        if not summary["tddft"][func]["covers_lam_min"]:
            print(f"WARNING: {a.nstates} states reach only {res['energies_ev'].max():.2f} eV; the MAC near "
                  f"{a.lam_min:.0f} nm is underestimated - increase --nstates.", flush=True)
        print(f"\n=== {func}: lowest 10 singlet states ===")
        print(lines_df.drop(columns="Transitions").head(10).to_string(index=False))
        print(f"MAC max in window: {summary['tddft'][func]['mac_max_m2_kg']:.1f} m^2/kg at "
              f"{summary['tddft'][func]['lambda_max_in_window_nm']:.0f} nm\n", flush=True)

    summary["total_time_s"] = round(time.time() - t_start, 1)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    print(f"Done in {summary['total_time_s']:.0f} s -> {outdir}")
    return summary


if __name__ == "__main__":
    main()
