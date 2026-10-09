"""
Physics emulator for the Sentinel-2 inversion.

A pixel's state is z = [log10 B, f_n, r, (dust)] where
    B   total algal abundance (cells mL^-1, BioSNICAR's unit, surface 2 cm layer)
    f_n fraction of Ancylonema nordenskioeldii cells (f_a = 1 - f_n is A. alaskanum),
        so the user-facing state [B, f_a, f_n, grain] lives on the simplex f_a + f_n = 1
    r   ice optical radius (granular grain radius, or bubble radius for solid ice), um
    dust optional mineral-dust nuisance (ppb, BioSNICAR's Greenland dust optics)

"Ours" (physics-informed): each species is a Phase 2 packaged cell of its own EMPIRICAL
geometry (Greenland volumes, Halbach et al. 2022; length:width, Prochazkova et al. 2021)
containing the phenolic pigment (Phase 1 TD-DFT Level 2 MAC calibrated against measured spectra of
the pigment - phase2/tddft_calibration.py - or the measured extract MAC of Williamson et al. 2020)
at the measured mass per cell scaled by the measured size dependence of concentration, plus
chlorophyll a/b and carotenoids (Williamson et al. 2020). Tier D adds the measured Fe-purpurogallin
absorption (Prochazkova et al. 2025) at the complexed fraction fitted to the S6 extract MAC.
Both species are mixed in BioSNICAR as two impurities.
"Tier A" (empirical baseline): BioSNICAR's default empirical glacier-algae optics,
which has no species information (f axis collapsed).

The emulator stores, at every grid node, the four Sentinel-2 band reflectances
(B2, B3, B4, B8; flux-weighted with the official ESA S2 spectral response functions), the
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
# bump when the forward physics changes, so cached emulators are rebuilt (v3: calibrated pigment MAC,
# Mie g, size-scaled species concentrations, measured clear-sky transmissivity; v4: meltwater-unit
# correction of cell counts, optional surface film)
PHYSICS_VERSION = "v4"
# dust nodes spanning the measured S6 loading (342 ug/g mean, 519 max; Cook et al. 2020), ppb
DUST_NODES_PPB = (3e4, 7e4, 1.5e5, 3e5, 6e5, 1.2e6)
S2_CENTRES_NM = (490, 560, 665, 842)


@dataclass
class SpeciesSpec:
    """Cell geometry and phenolic loading of one Ancylonema species (empirical defaults below)."""
    diameter_um: float
    length_um: float
    c_internal: float | None = None   # phenolic kg m^-3 of cell volume; None = empirical (Williamson 2020)


def empirical_species():
    """Greenland cell volumes (Halbach et al. 2022) with measured length:width ratios
    (Prochazkova et al. 2021); see phase2/empirical_data.species_geometry()."""
    import empirical_data as ED
    g = ED.species_geometry()
    return {k: SpeciesSpec(v["diameter_um"], v["length_um"]) for k, v in g.items()}


@dataclass
class EmulatorConfig:
    model: str = "ours"                       # 'ours' or 'tierA'
    tier: str = "C"                           # optical tier for 'ours': 'C' (uncomplexed) or 'D' (+ Fe complex)
    phenol: str = "tddft"                     # 'tddft' (calibrated Phase 1 Level 2), 'tddft_raw', 'williamson2020'
    photosynthetic: bool = True               # add chl a, chl b, carotenoids (Williamson et al. 2020)
    ice_mode: str = "bubbly"                  # empirical: field NIR requires solid bubbly ice (README)
    rho: float = 450.0                        # weathering crust (Cooper et al. 2018, mean 0.45 g cm-3)
    rho_bottom: float = 690.0                 # near-surface ice (Cooper et al. 2018, mean 0.69 g cm-3)
    sza: float = 47.0
    spacecraft: str = "S2A"                   # selects the ESA spectral response functions
    log_b: tuple = (1.0, 6.0, 0.25)           # start, stop, step (log10 cells/mL)
    f_n: tuple = (0.0, 1.0, 0.1)
    r_um: tuple = (300.0, 20000.0, 28)        # min, max, n (log-spaced; bubbly-ice optical radius)
    dust_ppb: tuple = ()                      # optional nuisance axis (ppb); () = no dust axis
    film_dz: float | None = None              # algae in a surface film of this thickness (m) inside the 2 cm
    film_only: bool = True                    # (see biosnicar_bridge.IceSpec); None = uniform over 2 cm
    sw_down: float | None = None              # broadband SW (W m^-2); None = clear-sky param.
    day_of_year: int = 196                    # for the Earth-Sun distance in the clear-sky SW
    species: dict = field(default_factory=empirical_species)

    def axes(self):
        ax = {"log_b": np.round(np.arange(self.log_b[0], self.log_b[1] + 1e-9, self.log_b[2]), 6)}
        ax["f_n"] = np.round(np.arange(self.f_n[0], self.f_n[1] + 1e-9, self.f_n[2]), 6) \
            if self.model == "ours" else np.array([0.5])
        ax["r_um"] = np.round(np.geomspace(self.r_um[0], self.r_um[1], int(self.r_um[2])), 3)
        if self.dust_ppb:
            ax["dust_ppb"] = np.array(self.dust_ppb, dtype=float)
        return ax


def _plain(x):
    """numpy scalars -> Python scalars so the metadata repr round-trips through ast.literal_eval."""
    if isinstance(x, dict):
        return {k: _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return type(x)(_plain(v) for v in x)
    return x.item() if isinstance(x, np.generic) else x


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
        """Atomic: written to a temporary file and renamed, so an interrupted save never leaves a
        truncated cache that a later run could load."""
        tmp = f"{path}.tmp{os.getpid()}.npz"
        np.savez_compressed(tmp, **{f"axis__{k}": v for k, v in self.axes.items()},
                            **{f"data__{k}": v for k, v in self.data.items()},
                            meta=np.array(repr(_plain(self.meta))))
        os.replace(tmp, path)

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
        import empirical_data as ED
        self.runner = bb.BioSNICARRunner(root, incoming=3)
        self.srf = ED.s2_srf_480(cfg.spacecraft, S2_BANDS)          # official ESA SRFs
        self.sw = cfg.sw_down if cfg.sw_down is not None else \
            bb.sw_down_clear_sky(cfg.sza, ED.clear_sky_transmissivity()[0], cfg.day_of_year)
        self.pg_per_cell = {}
        if cfg.model == "ours":
            kw = co.water_k_480(root)
            import tddft_calibration as TC
            if cfg.phenol == "williamson2020":
                if cfg.tier == "D":
                    raise ValueError("tier D is built from the calibrated TD-DFT spectrum (phenol='tddft')")
                ph = ED.pigment_macs_480()["phenolics_williamson2020"]
                mac = co.to_480(lambda wl: np.interp(wl, co.WVL_480_NM, ph))
            else:
                ligand = co.demo_spectrum(root) if demo else co.load_phase1(phase1_l2, "level2")
                if cfg.phenol == "tddft_raw":
                    if cfg.tier == "D":
                        raise ValueError("tier D needs the calibration (phenol='tddft')")
                    mac = co.to_480(ligand.mac_at)
                else:
                    cal = TC.cached_calibration(ligand, verbose=False)
                    mac = co.to_480(cal.mac_D if cfg.tier == "D" else cal.mac_C)
            self.imps = {}
            for name, sp in cfg.species.items():
                extra = co.empirical_pigments(species=name) if cfg.photosynthetic else ()
                c_i = sp.c_internal if sp.c_internal is not None else co.empirical_phenolic_concentration(name)
                cell = co.CellModel(CellGeometry("cylinder", sp.diameter_um / 2.0, sp.length_um),
                                    c_i, vd_diagnostic=False, extra_pigments=extra)
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
        return self.bb.IceSpec(r_um, self.cfg.rho, rho_bottom=self.cfg.rho_bottom, mode=self.cfg.ice_mode,
                               film_dz=self.cfg.film_dz, film_only=self.cfg.film_only)

    def _bands(self, alb, flx):
        """Flux-weighted band reflectance with the ESA spectral response functions."""
        out = []
        for b in S2_BANDS:
            w = self.srf[b] * flx
            out.append(float(np.sum(alb * w) / np.sum(w)))
        return out

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
    import biosnicar_bridge as bb
    import provenance as PV
    try:
        bsn = PV.biosnicar_revision(bb.locate_biosnicar(biosnicar))
    except ImportError:
        bsn = None
    # content fingerprint: configuration, Phase 1 result CONTENT (not its path), measured data, physics
    # code and BioSNICAR revision; PHYSICS_VERSION is kept as a human-readable label only
    tag = PV.fingerprint(dict(cfg={k: v for k, v in asdict(cfg).items() if k != "species"}, species=repr(cfg.species),
                              demo=demo, physics_version=PHYSICS_VERSION, phase1=PV.phase1_fingerprint(phase1_l2),
                              empirical=PV.empirical_data_fingerprint(), code=PV.code_fingerprint(),
                              biosnicar=bsn))
    if cache and os.path.isfile(cache):
        try:
            em = Emulator.load(cache)
        except (ValueError, SyntaxError, KeyError, OSError, EOFError, __import__("zipfile").BadZipFile):  # unreadable: rebuild
            em = None
        if em is not None and em.meta.get("tag") == tag:
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
