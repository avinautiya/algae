"""
Forward model for uncertainty propagation: one parameter vector -> albedo and
radiative forcing for tiers A-D. It is the Phase 2 physics (packaging, cell
optics, BioSNICAR bridge) with three speed-ups that do not change results:

  * Q* from a pre-computed look-up table (qstar_table.py, <0.3 % error);
  * clean-ice layer optics cached per BioSNICAR LUT radius (20 um steps);
  * the slow van Diedenhoven SSA diagnostic switched off (g is fixed at 0.96).

Outputs per sample
------------------
bba_clean, bba_{A,B,C,D}    broadband albedo 300-2500 nm
rf_{A,B,C,D}                instantaneous SW forcing vs clean ice [W m^-2]
d_CB = rf_C - rf_B          packaging effect (negative: packaging reduces forcing)
d_DC = rf_D - rf_C          Fe-complexation / aggregation effect
d_{B,C,D}A = rf_X - rf_A    anomaly relative to BioSNICAR's empirical algae
eff_{A,B,C,D}               forcing efficiency, rf / (conc / 10^4 cells mL^-1):
                            W m^-2 per 10^4 cells mL^-1 (removes the trivial
                            abundance scaling so micro-scale controls are visible)
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))
sys.path.insert(0, HERE)

TIERS = ("A", "B", "C", "D")


class ForwardModel:
    def __init__(self, phase1_l2: str | None = None, functional: str = "B3LYP", demo: bool = False,
                 biosnicar: str | None = None, sza: float = 45.0, sw_down: float | None = None,
                 transmissivity: float | None = None, incoming: int = 3, window=(350.0, 800.0),
                 uv_mode: str = "hold", ice_mode: str = "bubbly", rho_bottom: float = 690.0,
                 photosynthetic: bool = True,
                 dz_top: float = 0.02, dz_bottom: float = 2.0, g_fixed: float = 0.96,
                 qtable_cache: str | None = None, calibration_point: dict | None = None):
        import biosnicar_bridge as bb
        import cell_optics as co
        import empirical_data as ED
        import tddft_calibration as TC
        from qstar_table import QStarTable

        self.bb, self.co, self.ED, self.TC = bb, co, ED, TC
        self.root = bb.locate_biosnicar(biosnicar)
        self.ligand = co.demo_spectrum(self.root) if demo else co.load_phase1(phase1_l2, "level2", functional)
        self.demo = demo
        self.runner = bb.BioSNICARRunner(self.root, incoming=incoming)
        self.kw = co.water_k_480(self.root)
        self.qtab = QStarTable(cache=qtable_cache)
        self.tierA = self.runner.default_impurity("glacier_algae")
        self.sza = sza
        self.sw_fixed = sw_down
        self.tau = ED.clear_sky_transmissivity()[0] if transmissivity is None else transmissivity
        self.window, self.uv_mode = tuple(window), uv_mode
        self.rho_bottom = rho_bottom
        self.ice_kw = dict(rho_bottom=rho_bottom, dz_top=dz_top, dz_bottom=dz_bottom, mode=ice_mode)
        self.g_fixed = g_fixed
        # chlorophyll a/b + carotenoids at their measured per-cell concentrations (Williamson et al. 2020)
        self.extra = co.empirical_pigments() if photosynthetic else ()
        self.v_ref = ED.s6_biovolume_um3()[0]
        if calibration_point is None:
            calibration_point = TC.calibrate(self.ligand, verbose=False).point()
        self.cal = calibration_point          # defaults for absent molecular parameters

    def resolve(self, p: dict) -> dict:
        """Derived inputs of one sample: bubble radius from the ice SSA (at the bottom-layer density),
        size-scaled intracellular concentration, downwelling shortwave from the transmissivity."""
        q = dict(p)
        if "ice_ssa" in q:
            q["grain_um"] = float(self.ED.bubble_radius_um(q["ice_ssa"], self.rho_bottom))
        d, L = cell_dimensions(q)
        q["size_factor"] = (np.pi * d * d / 4.0 * L / self.v_ref) ** q.get("conc_size_exponent", 0.0)
        tau = q.get("transmissivity", self.tau)
        q["sw_down"] = self.sw_fixed if self.sw_fixed is not None else self.bb.sw_down_clear_sky(self.sza, tau)
        return q

    def warm(self, grain_min: float, grain_max: float):
        """Pre-compute clean-ice optics for every LUT radius in range (shared by forked workers)."""
        bb = self.bb
        if self.ice_kw["mode"] == "grains":   # granular optics depend on radius only -> cache them all
            for r in self.runner.available_radii("grains"):
                if grain_min - 20 <= r <= grain_max + 20:
                    self.runner.ice(bb.IceSpec(r, 650.0, **self.ice_kw))
        else:                                 # bubbly optics also depend on density: just load the tables
            self.runner.ice(bb.IceSpec(grain_min, 450.0, **self.ice_kw))
        self.runner.detach_luts()

    # ---- molecular -------------------------------------------------------------
    def mac480(self, p):
        """Calibrated pigment MAC (tier B/C) and with the Fe-complexed fraction (tier D); the
        molecular parameters are the calibration's dE (eV), FWHM (eV), f and phi."""
        c = self.cal
        dE, w = p.get("dE_ev", c["dE"]), p.get("fwhm_ev", c["w"])
        f, phi = p.get("f_scale", c["f"]), p.get("fe_fraction", c["phi"])
        mac_c = self.co.to_480(self.TC.perturbed_mac(self.ligand, dE, f, w), self.window, self.uv_mode)
        mac_d = self.co.to_480(self.TC.complexed_mac(self.ligand, dE, f, w, phi), self.window, self.uv_mode)
        return mac_c, mac_d

    # ---- cellular --------------------------------------------------------------
    def cell(self, p):
        from pigment_packaging import CellGeometry
        d, L = cell_dimensions(p)
        k = p.get("size_factor", 1.0)
        extra = tuple((m, c * k) for m, c in self.extra)
        return self.co.CellModel(CellGeometry("cylinder", d / 2.0, L), p["c_internal"] * k, q_func=self.qtab,
                                 vd_diagnostic=False, g_fixed=self.g_fixed, extra_pigments=extra)

    # ---- full chain ------------------------------------------------------------
    def evaluate(self, p: dict) -> dict:
        bb = self.bb
        p = self.resolve(p)
        mac_c, mac_d = self.mac480(p)
        cell = self.cell(p)
        oC = cell.optics(mac_c, self.kw, packaged=True)
        oB = cell.optics(mac_c, self.kw, packaged=False, scatter_from=oC)
        oD = cell.optics(mac_d, self.kw, packaged=True)
        imps = {"A": self.tierA,
                "B": bb.CustomImpurity("B", oB["ext_xsc"], oB["ss_alb"], oB["asm_prm"]),
                "C": bb.CustomImpurity("C", oC["ext_xsc"], oC["ss_alb"], oC["asm_prm"]),
                "D": bb.CustomImpurity("D", oD["ext_xsc"], oD["ss_alb"], oD["asm_prm"])}
        spec = bb.IceSpec(p["grain_um"], p["rho_top"], **self.ice_kw)
        conc = p["conc_cells_ml"]
        alb0, flx, _ = self.runner.run(spec, self.sza)
        out = {"bba_clean": self.runner.broadband(alb0, flx)}
        for t, imp in imps.items():
            alb, _, _ = self.runner.run(spec, self.sza, imp, conc)
            out[f"bba_{t}"] = self.runner.broadband(alb, flx)
            out[f"rf_{t}"] = self.runner.forcing(alb0, alb, flx, p["sw_down"])
            out[f"eff_{t}"] = out[f"rf_{t}"] / (conc / 1e4)
        out["d_CB"] = out["rf_C"] - out["rf_B"]
        out["d_DC"] = out["rf_D"] - out["rf_C"]
        for t in "BCD":
            out[f"d_{t}A"] = out[f"rf_{t}"] - out["rf_A"]
        out["pigment_pg_per_cell"] = cell.pigment_mass_per_cell_kg * 1e15
        out["cell_diameter_um"], out["cell_length_um"] = cell_dimensions(p)
        out["ice_radius_um"], out["sw_down"] = p["grain_um"], p["sw_down"]
        out["c_internal_eff"] = p["c_internal"] * p["size_factor"]
        return out


def cell_dimensions(p: dict):
    """(diameter, length) in um from either (cell_volume_um3, cell_aspect) or explicit dimensions."""
    if "cell_volume_um3" in p:
        ar = p["cell_aspect"]
        d = (4.0 * p["cell_volume_um3"] / (np.pi * ar)) ** (1.0 / 3.0)
        return d, ar * d
    return p["cell_diameter_um"], p["cell_length_um"]


# --------------------------------------------------------------------------- #
# Batch evaluation (optionally parallel)                                       #
# --------------------------------------------------------------------------- #
_WORKER = None


def _init_worker(kwargs):
    global _WORKER
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _WORKER = ForwardModel(**kwargs)


def _eval_chunk(rows):
    return [_WORKER.evaluate(r) for r in rows]


def evaluate_many(model_kwargs: dict, samples: list[dict], workers: int = 1, chunk: int = 64,
                  progress: bool = True):
    """Evaluate all parameter dicts; returns a pandas DataFrame (inputs + outputs).

    The model is built and warmed once in the parent; with the 'fork' start method
    (Linux / Colab) worker processes inherit it, including all caches.
    """
    import pandas as pd

    global _WORKER
    t0 = time.time()
    if _WORKER is None or getattr(_WORKER, "_kwargs", None) != model_kwargs:
        _WORKER = ForwardModel(**model_kwargs)
        _WORKER._kwargs = dict(model_kwargs)
    g = [_WORKER.resolve(s)["grain_um"] for s in samples]
    _WORKER.warm(min(g), max(g))
    chunks = [samples[i:i + chunk] for i in range(0, len(samples), chunk)]
    results = []
    import multiprocessing as mp
    if workers > 1 and "fork" in mp.get_all_start_methods():
        with mp.get_context("fork").Pool(workers) as pool:
            for k, res in enumerate(pool.imap(_eval_chunk, chunks)):
                results.extend(res)
                if progress:
                    _progress(k + 1, len(chunks), t0)
    elif workers > 1:     # spawn platforms: each worker rebuilds the model
        with mp.get_context().Pool(workers, initializer=_init_worker, initargs=(model_kwargs,)) as pool:
            for k, res in enumerate(pool.imap(_eval_chunk, chunks)):
                results.extend(res)
                if progress:
                    _progress(k + 1, len(chunks), t0)
    else:
        for k, c in enumerate(chunks):
            results.extend(_eval_chunk(c))
            if progress:
                _progress(k + 1, len(chunks), t0)
    df = pd.concat([pd.DataFrame(samples), pd.DataFrame(results)], axis=1)
    if progress:
        print(f"  {len(samples)} model runs in {time.time() - t0:.0f} s", flush=True)
    return df


def _progress(k, n, t0):
    if k == n or k % max(1, n // 10) == 0:
        el = time.time() - t0
        print(f"  {k}/{n} chunks  {el:6.0f} s elapsed, ~{el / k * (n - k):6.0f} s left", flush=True)
