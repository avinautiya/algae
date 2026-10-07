#!/usr/bin/env python
"""
Phase 2 driver: molecular MAC -> packaging-corrected cell optics -> BioSNICAR
albedo and radiative forcing for model tiers A-D -> Figures 2A-2C.

Examples
--------
# Normal run: reads Phase 1 results (results/level1, results/level2 next to phase1/)
python run_phase2.py --phase1-l1 ../phase1/results/level1 --phase1-l2 ../phase1/results/level2

# Pipeline check before Phase 1 is finished (BioSNICAR's ppg.csv as a stand-in; outputs marked DEMO)
python run_phase2.py --demo

# Tier D from real Level 3/4 TD-DFT output instead of the provisional surrogate
python run_phase2.py --level34-csv ../phase1/results/level3/level3_B3LYP_spectrum.csv --level34-molar-mass 908.5

Tiers
-----
A  BioSNICAR default empirical glacier-algae optics (ice_algae_empirical_Chevrollier2023)
B  Level 2 molecular MAC, pigment absorbing as if dissolved (no packaging)
C  Level 2 molecular MAC with Duysens packaging inside Ancylonema-like cells
D  Fe(III)-phenolic complexed + aggregated MAC (Level 3/4 CSV, or provisional surrogate), packaged
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

import biosnicar_bridge as bb  # noqa: E402


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    io = p.add_argument_group("inputs / outputs")
    io.add_argument("--phase1-l1", default=os.path.join(HERE, "..", "phase1", "results", "level1"))
    io.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    io.add_argument("--functional", default="B3LYP", help="which Phase 1 TD-DFT spectrum to use")
    io.add_argument("--demo", action="store_true", help="use BioSNICAR's ppg.csv instead of Phase 1 output")
    io.add_argument("--level34-csv", default=None, help="Level 3/4 MAC CSV (Wavelength_nm, MAC_estimated)")
    io.add_argument("--level34-molar-mass", type=float, default=None)
    io.add_argument("--biosnicar", default=None, help="path to a biosnicar-py checkout")
    io.add_argument("--outdir", default=os.path.join(HERE, "results"))
    io.add_argument("--usetex", action="store_true", help="render labels with a LaTeX installation")

    cell = p.add_argument_group("cell / packaging")
    cell.add_argument("--cell-shape", choices=["cylinder", "sphere"], default="cylinder")
    cell.add_argument("--cell-radius", type=float, default=None,
                      help="um (reference cell for BioSNICAR); default: A. nordenskioeldii, empirical_data")
    cell.add_argument("--cell-length", type=float, default=None, help="um (cylinders); default as above")
    cell.add_argument("--no-photosynthetic", action="store_true",
                      help="omit chlorophyll a/b and carotenoids (Williamson et al. 2020) from tiers B-D")
    cell.add_argument("--c-internal", type=float, default=None,
                      help="intracellular pigment concentration, kg m^-3 of cell volume")
    cell.add_argument("--c-internal-grid", type=float, nargs="+", default=None,
                      help="default: empirical c_i and c_i x (mean -/+ 1 SD of phenolics per cell)/mean")
    cell.add_argument("--vacuole-fraction", type=float, default=1.0)
    cell.add_argument("--sizes", type=float, nargs="+", default=None,
                      help="cell sizes (um) for the packaging grid / Fig. 2A")
    cell.add_argument("--size-is", choices=["length", "radius"], default="length",
                      help="interpret --sizes as cylinder length (diameter = length/aspect) or as radius")
    cell.add_argument("--aspect", type=float, default=2.0, help="cylinder length / diameter")
    cell.add_argument("--g-mode", choices=["fixed", "vd2014"], default="fixed")
    cell.add_argument("--window", type=float, nargs=2, default=[350.0, 800.0])
    cell.add_argument("--uv-mode", choices=["hold", "molecular", "zero"], default="hold")

    ice = p.add_argument_group("ice column / illumination sweep")
    ice.add_argument("--ice-mode", choices=["grains", "bubbly"], default="bubbly")
    ice.add_argument("--grains", type=float, nargs="+", default=[1000, 3000, 6000, 10000, 15000],
                     help="grain (or bubble) effective radius, um")
    ice.add_argument("--densities", type=float, nargs="+", default=[330, 450, 560],
                     help="surface-layer density, kg m^-3")
    ice.add_argument("--rho-bottom", type=float, default=690.0,
                      help="Cooper et al. (2018) near-surface ice mean, kg/m3")
    ice.add_argument("--dz-top", type=float, default=0.02, help="algae-bearing layer thickness, m")
    ice.add_argument("--dz-bottom", type=float, default=2.0)
    ice.add_argument("--sza", type=float, nargs="+", default=[45, 55, 65, 75])
    ice.add_argument("--conc-min", type=float, default=1e3)
    ice.add_argument("--conc-max", type=float, default=1e7)
    ice.add_argument("--n-conc", type=int, default=17)
    ice.add_argument("--incoming", type=int, default=3, help="BioSNICAR spectrum: 3 = sub-Arctic summer")
    ice.add_argument("--diffuse", action="store_true", help="diffuse (cloudy) instead of direct beam")
    ice.add_argument("--sw-down", type=float, default=None,
                     help="fixed broadband SW down (W m^-2); default: clear-sky parameterisation per SZA")
    ice.add_argument("--transmissivity", type=float, default=0.75)

    ref = p.add_argument_group("reference state for figures")
    ref.add_argument("--ref-grain", type=float, default=10000,
                     help="um; Cooper et al. (2021) bubbly-ice optical radius 9.3-10.6 mm at S6")
    ref.add_argument("--ref-density", type=float, default=450,
                     help="kg/m3; Cooper et al. (2018) weathering-crust mean")
    ref.add_argument("--ref-sza", type=float, default=45, help="deg; ~solar noon at S6 in July")
    ref.add_argument("--fig2c-concs", type=float, nargs="+", default=[1e4, 1e5])
    ref.add_argument("--baseline", default="A", choices=list("ABCD"), help="tier subtracted in Fig. 2C")
    return p.parse_args(argv)


def _empirical_defaults(a):
    """Fill unset cell parameters from published measurements (phase2/empirical_data.py)."""
    import empirical_data as ED
    g = ED.species_geometry()["nordenskioeldii"]
    if a.cell_radius is None:
        a.cell_radius = g["diameter_um"] / 2.0
    if a.cell_length is None:
        a.cell_length = g["length_um"]
    ci = ED.intracellular_concentration_kg_m3("phenolics")
    if a.c_internal is None:
        a.c_internal = ci
    if a.c_internal_grid is None:
        m, sd = ED.pigments_per_cell()["phenolics"]
        a.c_internal_grid = [ci * (m - sd) / m, ci, ci * (m + sd) / m]


def size_to_geom(size_um, a, CellGeometry):
    if a.cell_shape == "sphere":
        return CellGeometry("sphere", size_um if a.size_is == "radius" else size_um / 2.0)
    if a.size_is == "length":
        return CellGeometry("cylinder", size_um / a.aspect / 2.0, size_um)
    return CellGeometry("cylinder", size_um, 2.0 * size_um * a.aspect)


def main(argv=None):
    a = parse_args(argv)
    t0 = time.time()
    root = bb.locate_biosnicar(a.biosnicar)
    _empirical_defaults(a)

    import cell_optics as co
    import figures as F
    from pigment_packaging import CellGeometry, mac_vivo

    tabdir, figdir = os.path.join(a.outdir, "tables"), os.path.join(a.outdir, "figures")
    os.makedirs(tabdir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)
    F.set_style(usetex=a.usetex)

    # ------------------------------------------------------------------ inputs
    if a.demo:
        l2 = co.demo_spectrum(root)
        l1 = None
        print("DEMO MODE: molecular MAC = BioSNICAR data/pigments/ppg.csv (not Phase 1 output)")
    else:
        l2 = co.load_phase1(a.phase1_l2, "level2", a.functional)
        try:
            l1 = co.load_phase1(a.phase1_l1, "level1", a.functional)
        except FileNotFoundError:
            l1 = None
            print("Level 1 results not found - Fig. 2A shows Level 2 only")
    print(f"Level 2 pigment: {l2.name} [{l2.source}]")

    if a.level34_csv:
        mm = a.level34_molar_mass or l2.molar_mass
        l34 = co.load_mac_csv(a.level34_csv, "Level 3/4", mm)
        mac_D_fn = l34.mac_at
        d_source = f"Level 3/4 CSV {a.level34_csv}"
        surrogate = None
    else:
        surrogate = co.FePhenolicSurrogate()
        mac_D_fn = lambda wl: surrogate.mac(l2, wl)  # noqa: E731
        d_source = "PROVISIONAL Fe(III)-phenolic surrogate (replace with Level 3/4 TD-DFT)"
    print(f"Tier D source: {d_source}")

    # --------------------------------------------- Task 1: packaging (300-800 nm)
    wl = np.arange(300.0, 801.0, 1.0)
    mac_l2 = l2.mac_at(wl)
    if a.sizes is None:      # the two Ancylonema species (Chevrollier 2022 volume, Prochazkova 2021 shape)
        import empirical_data as ED
        sizes_geom = [CellGeometry("cylinder", g["diameter_um"] / 2.0, g["length_um"])
                      for g in ED.species_geometry().values()]
    else:
        sizes_geom = [size_to_geom(s, a, CellGeometry) for s in a.sizes]
    grid = np.empty((len(sizes_geom), len(a.c_internal_grid), wl.size))
    qgrid = np.empty_like(grid)
    for i, g in enumerate(sizes_geom):
        for j, c in enumerate(a.c_internal_grid):
            grid[i, j], qgrid[i, j], _ = mac_vivo(mac_l2, c, g, a.vacuole_fraction)
    np.savez(os.path.join(tabdir, "mac_vivo_grid.npz"), wavelength_nm=wl, mac_raw=mac_l2,
             mac_vivo=grid, q_star=qgrid, c_internal_kg_m3=np.array(a.c_internal_grid),
             cell_radius_um=np.array([g.radius for g in sizes_geom]),
             cell_length_um=np.array([g.length for g in sizes_geom]), cell_shape=a.cell_shape)
    long = [dict(cell=g.label(plain=True), radius_um=g.radius, length_um=g.length,
                 c_internal_kg_m3=c, Wavelength_nm=w, MAC_raw=mr, MAC_vivo=mv, Q_star=q)
            for i, g in enumerate(sizes_geom) for j, c in enumerate(a.c_internal_grid)
            for w, mr, mv, q in zip(wl, mac_l2, grid[i, j], qgrid[i, j])]
    pd.DataFrame(long).to_csv(os.path.join(tabdir, "mac_vivo_grid.csv"), index=False, float_format="%.6g")

    # Fig 2A at the reference concentration
    jc = int(np.argmin(np.abs(np.array(a.c_internal_grid) - a.c_internal)))
    if not np.isclose(a.c_internal_grid[jc], a.c_internal):
        mv_ref = {g.label(): mac_vivo(mac_l2, a.c_internal, g, a.vacuole_fraction)[0] for g in sizes_geom}
        q_ref = {g.label(): mac_vivo(mac_l2, a.c_internal, g, a.vacuole_fraction)[1] for g in sizes_geom}
    else:
        mv_ref = {g.label(): grid[i, jc] for i, g in enumerate(sizes_geom)}
        q_ref = {g.label(): qgrid[i, jc] for i, g in enumerate(sizes_geom)}
    raw = {"Level 2 glucoside": mac_l2} if l1 is None else {"Level 2 glucoside": mac_l2, "Level 1 core": l1.mac_at(wl)}
    pck = np.loadtxt(os.path.join(root, "data", "pigments", "pckg_GA.csv"))
    wl_p = 200.0 + np.arange(pck.size)
    keep = (wl_p >= 300) & (wl_p <= 600)      # empirical factor exceeds 1 (noise) beyond ~600 nm
    fig = F.fig2a(wl, raw, mv_ref, q_ref, figdir, empirical_q=(wl_p[keep], pck[keep]),
                  c_internal=a.c_internal, demo=a.demo)
    F.save(fig, figdir, "Fig2A_MAC_packaging")

    # ------------------------------------------------ cell optics, 480 bands
    kw = co.water_k_480(root)
    ref_geom = CellGeometry(a.cell_shape, a.cell_radius, a.cell_length)
    extra = () if a.no_photosynthetic else co.empirical_pigments()
    cell = co.CellModel(ref_geom, a.c_internal, a.vacuole_fraction, g_mode=a.g_mode, extra_pigments=extra)
    mac480_L2 = co.to_480(l2.mac_at, tuple(a.window), a.uv_mode)
    mac480_D = co.to_480(mac_D_fn, tuple(a.window), a.uv_mode)
    optics = {"A": co.model_a_optics(root)}
    optics["C"] = cell.optics(mac480_L2, kw, packaged=True)
    optics["B"] = cell.optics(mac480_L2, kw, packaged=False, scatter_from=optics["C"])
    optics["D"] = cell.optics(mac480_D, kw, packaged=True)
    optics = {k: optics[k] for k in "ABCD"}

    rows = {"Wavelength_nm": co.WVL_480_NM}
    for t, o in optics.items():
        for k in ("ext_xsc", "abs_xsc", "ss_alb", "asm_prm"):
            rows[f"{t}_{k}"] = o[k]
    rows["MAC_L2_480"], rows["MAC_D_480"] = mac480_L2, mac480_D
    rows["C_q_star"] = optics["C"]["q_star"]
    rows["C_ssa_vandiedenhoven2014"] = optics["C"]["ssa_vandiedenhoven"]
    pd.DataFrame(rows).to_csv(os.path.join(tabdir, "cell_optics_480band.csv"), index=False, float_format="%.6g")
    os.makedirs(os.path.join(a.outdir, "lap_entries"), exist_ok=True)
    for t in "BCD":
        bb.export_lap_entry(os.path.join(a.outdir, "lap_entries", f"tier{t}_glacier_algae.npz"),
                            f"tier{t}_glacier_algae", optics[t])
    F.save(F.fig_s1_cell_optics(co.WVL_480_NM, optics, figdir), figdir, "FigS1_cell_optics")

    vis = (co.WVL_480_NM >= 350) & (co.WVL_480_NM <= 700)
    dssa = np.abs(optics["C"]["ss_alb"][vis] - optics["C"]["ssa_vandiedenhoven"][vis]).max()
    print(f"Reference cell: {ref_geom.label(plain=True)}, c_i = {a.c_internal:g} kg/m3 -> "
          f"{cell.pigment_mass_per_cell_kg * 1e15:.1f} pg pigment/cell; "
          f"max |SSA(Duysens) - SSA(van Diedenhoven 2014)| (350-700 nm) = {dssa:.3f}")

    # ------------------------------------------------ Task 2/3: BioSNICAR sweep
    runner = bb.BioSNICARRunner(root, incoming=a.incoming, direct=0 if a.diffuse else 1)
    imps = {"A": runner.default_impurity("glacier_algae")}
    for t in "BCD":
        imps[t] = bb.CustomImpurity(f"tier{t}", optics[t]["ext_xsc"], optics[t]["ss_alb"], optics[t]["asm_prm"])
    concs = np.logspace(np.log10(a.conc_min), np.log10(a.conc_max), a.n_conc)
    concs = np.unique(np.concatenate([concs, a.fig2c_concs]))
    records, spectra_ref = [], {}
    for grain in a.grains:
        for rho in a.densities:
            spec = bb.IceSpec(grain, rho, a.rho_bottom, a.dz_top, a.dz_bottom, a.ice_mode)
            for sza in a.sza:
                sw = a.sw_down if a.sw_down is not None else bb.sw_down_clear_sky(sza, a.transmissivity)
                alb0, flx, bba0 = runner.run(spec, sza)
                base = dict(grain_um=grain, rho_top=rho, sza=sza, sw_down=sw)
                records.append(dict(base, tier="clean", conc=0.0, bba=runner.broadband(alb0, flx),
                                    bba_biosnicar=bba0, rf=0.0))
                is_ref = (grain == a.ref_grain and rho == a.ref_density and sza == a.ref_sza)
                if is_ref:
                    spectra_ref["clean"] = alb0
                for t, imp in imps.items():
                    for c in concs:
                        alb, _, bba_full = runner.run(spec, sza, imp, c)
                        records.append(dict(base, tier=t, conc=c, bba=runner.broadband(alb, flx),
                                            bba_biosnicar=bba_full, rf=runner.forcing(alb0, alb, flx, sw)))
                        if is_ref and np.isclose(c, a.fig2c_concs[0]):
                            spectra_ref[t] = alb
    df = pd.DataFrame.from_records(records)
    df.to_csv(os.path.join(tabdir, "biosnicar_sweep.csv"), index=False, float_format="%.6g")
    np.savez(os.path.join(tabdir, "spectral_albedo_reference.npz"), wavelength_nm=co.WVL_480_NM,
             **{k: v for k, v in spectra_ref.items()})
    print(f"BioSNICAR sweep: {len(df)} runs in {time.time() - t0:.0f} s total elapsed")

    # ------------------------------------------------ Task 3 summary table
    ref = dict(grain_um=a.ref_grain, rho_top=a.ref_density, sza=a.ref_sza, ice_mode=a.ice_mode)
    if not (a.ref_grain in a.grains and a.ref_density in a.densities and a.ref_sza in a.sza):
        raise SystemExit("reference state must be one of the swept grains/densities/SZAs")
    r = df[(df.grain_um == a.ref_grain) & (df.rho_top == a.ref_density) & (df.sza == a.ref_sza)]
    summ = (r[r.conc.isin(a.fig2c_concs) | (r.tier == "clean")]
            .pivot_table(index="tier", columns="conc", values=["bba", "rf"]).round(4))
    summ.to_csv(os.path.join(tabdir, "summary_reference_state.csv"))
    print("\nReference state: radius %.0f um, rho %.0f kg/m3, SZA %.0f deg, SW_down %.0f W/m2"
          % (a.ref_grain, a.ref_density, a.ref_sza, r.sw_down.iloc[0]))
    print(summ.to_string())

    # ------------------------------------------------ Task 4: figures
    F.save(F.fig2b(df, ref, figdir, demo=a.demo), figdir, "Fig2B_albedo_vs_concentration")
    F.save(F.fig2c(df, ref, a.fig2c_concs, figdir, baseline=a.baseline, demo=a.demo), figdir,
           "Fig2C_forcing_anomaly_vs_grain")
    if spectra_ref:
        F.save(F.fig_s2_spectral_albedo(co.WVL_480_NM, {k: v for k, v in spectra_ref.items() if k != "clean"},
                                        spectra_ref["clean"], a.fig2c_concs[0], figdir),
               figdir, "FigS2_spectral_albedo")

    cfg = vars(a).copy()
    cfg.update(demo=a.demo, level2_source=l2.source, tier_d_source=d_source,
               tier_d_surrogate=(surrogate.__dict__ if surrogate else None),
               reference_cell=co.describe(cell), biosnicar_root=root,
               runtime_s=round(time.time() - t0, 1))
    with open(os.path.join(a.outdir, "run_config.json"), "w") as fh:
        json.dump(cfg, fh, indent=2, default=str)
    print(f"\nDone in {time.time() - t0:.0f} s -> {a.outdir}")
    return df


if __name__ == "__main__":
    main()
