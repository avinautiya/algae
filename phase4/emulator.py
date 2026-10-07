"""
Physics emulator for the Sentinel-2 inversion.

A pixel's state is z = [log10 B, f_n, r, (dust)] where
    B   total algal abundance (cells mL^-1, BioSNICAR's unit, surface 2 cm layer)
    f_n fraction of Ancylonema nordenskioeldii cells (f_a = 1 - f_n is A. alaskanum),
        so the user-facing state [B, f_a, f_n, grain] lives on the simplex f_a + f_n = 1
    r   ice optical radius (granular grain radius, or bubble radius for solid ice), um
    dust optional mineral-dust nuisance (ppb, BioSNICAR's Greenland dust optics)

"Ours" (physics-informed): each species is a Phase 2 packaged cell (tier C, or tier D
with the Fe-phenolic surrogate) of its own geometry, built from the Phase 1 TD-DFT
spectrum. Both species are mixed in BioSNICAR as two impurities.
"Tier A" (empirical baseline): BioSNICAR's default empirical glacier-algae optics,
which has no species information (f axis collapsed).

The emulator stores, at every grid node, the four Sentinel-2 band reflectances
(B2, B3, B4, B8; flux-weighted SRF convolution with BioSNICAR's S2 SRFs), the
broadband albedo (300-2500 nm), the instantaneous forcing of algae and of all
impurities, and the pigment mass per mL. It is then refined along log B with
monotone cubic interpolation for inference.
"""

from __future__ import annotations

import copy
import os
import sys
import time
from dataclasses import dataclass, field, asdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

S2_BANDS = ("B2", "B3", "B4", "B8")
S2_CENTRES_NM = (490, 560, 665, 842)


@dataclass
class SpeciesSpec:
    """Cell geometry and pigment loading of one Ancylonema species. These defaults are
    placeholders within published size ranges - set them from your own microscopy."""
    diameter_um: float
    length_um: float
    c_internal: float = 50.0          # kg m^-3 of cell volume (Phase 2 default)


DEFAULT_SPECIES = {
    "nordenskioeldii": SpeciesSpec(diameter_um=10.0, length_um=22.0),   # cells of filaments
    "alaskanum": SpeciesSpec(diameter_um=11.0, length_um=16.0),          # shorter single cells
}


@dataclass
class EmulatorConfig:
    model: str = "ours"                       # 'ours' or 'tierA'
    tier: str = "C"                           # optical tier for 'ours': 'C' or 'D'
    ice_mode: str = "grains"                  # 'grains' (granular) or 'bubbly' (solid ice)
    rho: float = 650.0                        # surface-layer density, kg m^-3
    sza: float = 47.0
    log_b: tuple = (1.0, 6.0, 0.2)            # start, stop, step (log10 cells/mL)
    f_n: tuple = (0.0, 1.0, 0.1)
    r_um: tuple = (1000.0, 3000.0, 80.0)
    dust_ppb: tuple = ()                      # e.g. (0, 1e3, 3e3, 1e4, 3e4, 1e5); () = no dust axis
    sw_down: float | None = None              # broadband SW (W m^-2); None = clear-sky param.
    species: dict = field(default_factory=lambda: dict(DEFAULT_SPECIES))

    def axes(self):
        ax = {"log_b": np.round(np.arange(self.log_b[0], self.log_b[1] + 1e-9, self.log_b[2]), 6)}
        ax["f_n"] = np.round(np.arange(self.f_n[0], self.f_n[1] + 1e-9, self.f_n[2]), 6) \
            if self.model == "ours" else np.array([0.5])
        ax["r_um"] = np.round(np.arange(self.r_um[0], self.r_um[1] + 1e-9, self.r_um[2]), 3)
        if self.dust_ppb:
            ax["dust_ppb"] = np.array(self.dust_ppb, dtype=float)
        return ax


class Emulator:
    """Grid of forward-model outputs; axes in a fixed order (log_b, f_n, r_um[, dust_ppb])."""

    OUTPUTS = ("bands", "bba", "rf_algae", "rf_total", "pigment_ug_l")

    def __init__(self, axes: dict, data: dict, meta: dict):
        self.axes = {k: np.asarray(v, dtype=float) for k, v in axes.items()}
        self.data = data
        self.meta = meta

    @property
    def names(self):
        return list(self.axes)

    @property
    def shape(self):
        return tuple(len(v) for v in self.axes.values())

    # ---- persistence -------------------------------------------------------------
    def save(self, path):
        np.savez_compressed(path, **{f"axis__{k}": v for k, v in self.axes.items()},
                            **{f"data__{k}": v for k, v in self.data.items()},
                            meta=np.array(repr(self.meta)))

    @classmethod
    def load(cls, path):
        d = np.load(path, allow_pickle=False)
        axes = {k[6:]: d[k] for k in d.files if k.startswith("axis__")}
        order = [n for n in ("log_b", "f_n", "r_um", "dust_ppb") if n in axes]
        data = {k[6:]: d[k] for k in d.files if k.startswith("data__")}
        import ast
        return cls({n: axes[n] for n in order}, data, ast.literal_eval(str(d["meta"])))

    # ---- refinement --------------------------------------------------------------
    def refine_log_b(self, step: float = 0.05):
        """Monotone cubic (PCHIP) interpolation along log B onto a finer grid."""
        from scipy.interpolate import PchipInterpolator
        lb = self.axes["log_b"]
        new = np.round(np.arange(lb[0], lb[-1] + 1e-9, step), 6)
        data = {k: (PchipInterpolator(lb, v, axis=0)(new) if np.all(np.isfinite(v))
                    else np.full((len(new),) + v.shape[1:], np.nan)) for k, v in self.data.items()}
        axes = dict(self.axes)
        axes["log_b"] = new
        return Emulator(axes, data, dict(self.meta, refined_step=step))

    @property
    def active(self):
        """Axes with more than one node (Tier A has a single, inactive f_n node)."""
        return [n for n in self.names if len(self.axes[n]) > 1]

    def coord(self, name):
        """Interpolation coordinate of an axis (dust is interpolated in log10(dust + 100))."""
        v = self.axes[name]
        return np.log10(v + 100.0) if name == "dust_ppb" else v

    def interpolator(self, key: str = "bands"):
        """Continuous (multilinear) interpolator over the ACTIVE axes, for MCMC."""
        from scipy.interpolate import RegularGridInterpolator
        sq = tuple(i for i, n in enumerate(self.names) if len(self.axes[n]) == 1)
        vals = np.squeeze(self.data[key], axis=sq) if sq else self.data[key]
        return RegularGridInterpolator(tuple(self.coord(n) for n in self.active), vals,
                                       bounds_error=False, fill_value=None)


# --------------------------------------------------------------------------- #
# Building                                                                      #
# --------------------------------------------------------------------------- #
class _Builder:
    """Holds BioSNICAR objects and species optics; evaluates one node."""

    def __init__(self, cfg: EmulatorConfig, phase1_l2: str | None, demo: bool, biosnicar: str | None):
        import biosnicar_bridge as bb
        import cell_optics as co
        from pigment_packaging import CellGeometry

        self.bb, self.cfg = bb, cfg
        root = bb.locate_biosnicar(biosnicar)
        from biosnicar.bands._core import load_srf
        self.runner = bb.BioSNICARRunner(root, incoming=3)
        self.srf = load_srf("sentinel2_msi")
        self.sw = cfg.sw_down if cfg.sw_down is not None else bb.sw_down_clear_sky(cfg.sza)
        self.pg_per_cell = {}
        if cfg.model == "ours":
            ligand = co.demo_spectrum(root) if demo else co.load_phase1(phase1_l2, "level2")
            kw = co.water_k_480(root)
            if cfg.tier == "D":
                sur = co.FePhenolicSurrogate()
                mac = co.to_480(lambda wl: sur.mac(ligand, wl))
            else:
                mac = co.to_480(ligand.mac_at)
            self.imps = {}
            for name, sp in cfg.species.items():
                cell = co.CellModel(CellGeometry("cylinder", sp.diameter_um / 2.0, sp.length_um),
                                    sp.c_internal, vd_diagnostic=False)
                o = cell.optics(mac, kw, packaged=True)
                self.imps[name] = bb.CustomImpurity(name, o["ext_xsc"], o["ss_alb"], o["asm_prm"])
                self.pg_per_cell[name] = cell.pigment_mass_per_cell_kg * 1e15
        else:
            self.imps = {"tierA": self.runner.default_impurity("glacier_algae")}
        d = np.load(os.path.join(root, "data", "OP_data", "480band", "lap.npz"))
        st = "dust_greenland_Cook_CENTRAL_20190911"
        self.dust = copy.deepcopy(self.runner.default_impurity("glacier_algae"))
        self.dust.name, self.dust.unit = "dust", 0
        self.dust.mac, self.dust.ssa, self.dust.g = d[st + "__ext_cff_mss"], d[st + "__ss_alb"], d[st + "__asm_prm"]
        self._clean = {}

    def spec(self, r_um):
        return self.bb.IceSpec(r_um, self.cfg.rho, mode=self.cfg.ice_mode)

    def _bands(self, alb, flx):
        from biosnicar.bands._core import srf_convolve
        return [srf_convolve(alb, flx, self.srf[b]) for b in S2_BANDS]

    def node(self, log_b, f_n, r_um, dust):
        B = 10.0 ** log_b
        spec = self.spec(r_um)
        if self.cfg.model == "ours":
            imps = [(self.imps["nordenskioeldii"], B * f_n), (self.imps["alaskanum"], B * (1.0 - f_n))]
            pig = B * (f_n * self.pg_per_cell["nordenskioeldii"] + (1 - f_n) * self.pg_per_cell["alaskanum"]) * 1e-3
        else:
            imps = [(self.imps["tierA"], B)]
            pig = np.nan
        key = (r_um, dust)
        if key not in self._clean:
            a0, flx, _ = self.runner.run_multi(spec, self.cfg.sza, [])
            ad, _, _ = self.runner.run_multi(spec, self.cfg.sza, [(self.dust, dust)]) if dust > 0 else (a0, flx, 0)
            self._clean[key] = (a0, ad, flx)
        a0, ad, flx = self._clean[key]
        alb, _, _ = self.runner.run_multi(spec, self.cfg.sza, imps + ([(self.dust, dust)] if dust > 0 else []))
        return (self._bands(alb, flx), self.runner.broadband(alb, flx),
                self.runner.forcing(ad, alb, flx, self.sw), self.runner.forcing(a0, alb, flx, self.sw), pig)


_BUILDER = None


def _build_r(args):
    ir, r_um, axes = args
    b = _BUILDER
    lb, fn = axes["log_b"], axes["f_n"]
    du = axes.get("dust_ppb", np.array([0.0]))
    out = np.empty((len(lb), len(fn), len(du), 8))
    for i, x in enumerate(lb):
        for j, f in enumerate(fn):
            for k, d in enumerate(du):
                bands, bba, rfa, rft, pig = b.node(x, f, r_um, d)
                out[i, j, k] = [*bands, bba, rfa, rft, pig]
    return ir, out


def build_emulator(cfg: EmulatorConfig, phase1_l2: str | None = None, demo: bool = False,
                   biosnicar: str | None = None, workers: int = 1, cache: str | None = None,
                   verbose: bool = True) -> Emulator:
    global _BUILDER
    tag = repr({k: v for k, v in asdict(cfg).items() if k != "species"}) + repr(cfg.species) + str(demo)
    if cache and os.path.isfile(cache):
        em = Emulator.load(cache)
        if em.meta.get("tag") == tag:
            if verbose:
                print(f"Emulator ({cfg.model}) loaded from {cache}")
            return em
    t0 = time.time()
    _BUILDER = _Builder(cfg, phase1_l2, demo, biosnicar)
    axes = cfg.axes()
    # snap r to BioSNICAR's look-up-table radii so nodes are exact model states
    axes["r_um"] = np.unique([_BUILDER.runner.snap_radius(r, cfg.ice_mode) for r in axes["r_um"]]).astype(float)
    nb, nf, nr = len(axes["log_b"]), len(axes["f_n"]), len(axes["r_um"])
    nd = len(axes.get("dust_ppb", [0]))
    if verbose:
        print(f"Building emulator ({cfg.model}, tier {cfg.tier if cfg.model == 'ours' else 'A'}, "
              f"{cfg.ice_mode}): {nb}x{nf}x{nr}x{nd} = {nb * nf * nr * nd} BioSNICAR runs", flush=True)
    raw = np.empty((nb, nf, nd, nr, 8))
    jobs = [(ir, r, axes) for ir, r in enumerate(axes["r_um"])]
    # Do ALL file I/O in the parent before forking: BioSNICAR keeps lazily-read .npz
    # handles at module level, and forked workers sharing one file offset corrupt reads.
    _BUILDER.runner.illumination(cfg.sza)
    for r in axes["r_um"]:
        _BUILDER.runner.ice(_BUILDER.spec(r))
    import multiprocessing as mp
    if workers > 1 and "fork" in mp.get_all_start_methods():
        with mp.get_context("fork").Pool(workers) as pool:
            for ir, out in pool.imap_unordered(_build_r, jobs):
                raw[:, :, :, ir] = out
    else:
        for job in jobs:
            ir, out = _build_r(job)
            raw[:, :, :, ir] = out
    raw = np.moveaxis(raw, 2, 3)                     # (b, f, r, d, 8)
    if "dust_ppb" not in axes:
        raw = raw[:, :, :, 0]
    data = {"bands": raw[..., :4], "bba": raw[..., 4], "rf_algae": raw[..., 5],
            "rf_total": raw[..., 6], "pigment_ug_l": raw[..., 7]}
    meta = dict(tag=tag, model=cfg.model, tier=cfg.tier, ice_mode=cfg.ice_mode, rho=cfg.rho, sza=cfg.sza,
                sw_down=_BUILDER.sw, pg_per_cell=_BUILDER.pg_per_cell, build_s=round(time.time() - t0, 1))
    em = Emulator(axes, data, meta)
    if cache:
        em.save(cache)
    if verbose:
        print(f"  built in {time.time() - t0:.0f} s", flush=True)
    return em


def direct_forward_chunk(states: list[dict]):
    """direct_forward with the already-initialised global builder (for forked workers)."""
    return direct_forward(_BUILDER.cfg, states)


def direct_forward(cfg: EmulatorConfig, states: list[dict], phase1_l2=None, demo=False, biosnicar=None):
    """Exact BioSNICAR evaluation at arbitrary (off-grid) states - used to synthesise
    observations and to test the emulator. Returns array (n, 4) and dict of extras."""
    global _BUILDER
    if _BUILDER is None or _BUILDER.cfg is not cfg:
        _BUILDER = _Builder(cfg, phase1_l2, demo, biosnicar)
    out = np.empty((len(states), 8))
    for i, s in enumerate(states):
        r = _BUILDER.runner.snap_radius(s["r_um"], cfg.ice_mode)
        bands, bba, rfa, rft, pig = _BUILDER.node(s["log_b"], s.get("f_n", 0.5), r, s.get("dust_ppb", 0.0))
        out[i] = [*bands, bba, rfa, rft, pig]
    return out
