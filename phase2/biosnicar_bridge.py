"""
Bridge between our cell optics and the BioSNICAR radiative-transfer model.

BioSNICAR (github.com/jmcook1186/biosnicar-py, v2.x) resolves its data folder
relative to the repository checkout, so it is used from a git clone placed on
sys.path (see locate_biosnicar / the Colab notebook), not from a pip wheel.

Integration strategy
--------------------
We do NOT edit BioSNICAR's data archive. Its impurity mixing routine
(`column_OPs.mix_in_impurities`) only needs objects exposing
`.unit, .conc, .mac, .ssa, .g`; we hand it `CustomImpurity` instances whose
arrays are our 480-band per-cell optics, then call BioSNICAR's own
adding-doubling solver. Model A uses BioSNICAR's own Impurity class and file,
so all four tiers pass through identical code.

`export_lap_entry` additionally writes our arrays in BioSNICAR's lap.npz key
convention (`<stem>__ext_xsc`, `__ss_alb`, `__asm_prm`) if you want to register
them permanently (see biosnicar-py/ADDING_DATA.md).
"""

from __future__ import annotations

import copy
import os
import sys
from dataclasses import dataclass

import numpy as np

S0 = 1361.0                     # total solar irradiance, W m^-2 (Kopp & Lean 2011)


def locate_biosnicar(path: str | None = None) -> str:
    """Put a BioSNICAR checkout on sys.path and return its root."""
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [path, os.environ.get("BIOSNICAR_PATH"),
                  os.path.join(here, "biosnicar-py"), os.path.join(here, "..", "biosnicar-py"),
                  os.path.join(here, "..", "..", "biosnicar-py"), "/content/biosnicar-py"]
    for c in candidates:
        if c and os.path.isfile(os.path.join(c, "biosnicar", "inputs.yaml")):
            root = os.path.abspath(c)
            if root not in sys.path:
                sys.path.insert(0, root)
            import biosnicar  # noqa: F401
            return root
    raise ImportError(
        "BioSNICAR checkout not found. Run:\n"
        "  git clone --depth 1 https://github.com/jmcook1186/biosnicar-py.git\n"
        "next to this repository (or set BIOSNICAR_PATH / --biosnicar).")


class CustomImpurity:
    """Duck-typed stand-in for biosnicar.classes.impurity.Impurity (cells/mL)."""

    def __init__(self, name: str, ext_xsc, ss_alb, asm_prm, conc=None):
        self.name = name
        self.file = f"{name} (in-memory)"
        self.unit = 1                   # 1 = cells per mL -> mac is m^2 per cell
        self.conc = conc or [0.0]
        self.mac = np.asarray(ext_xsc, dtype=float)
        self.ssa = np.asarray(ss_alb, dtype=float)
        self.g = np.asarray(asm_prm, dtype=float)
        assert self.mac.shape == self.ssa.shape == self.g.shape == (480,)


@dataclass
class IceSpec:
    """Bare-ice column: algae-bearing weathering crust (dz_top) over clean ice (dz_bottom); with film_dz
    the crust is split into two sub-layers (3 layers in total).

    mode='grains': granular ice spheres of radius `rds_um` (layer_type 0)
    mode='bubbly': solid ice with Fresnel surface, `rds_um` = bubble radius (layer_type 1)
    """
    rds_um: float = 1500.0
    rho_top: float = 650.0
    rho_bottom: float = 850.0
    dz_top: float = 0.02
    dz_bottom: float = 2.0
    mode: str = "grains"
    film_dz: float | None = None       # if set (< dz_top), the top layer is split into [film_dz, dz_top - film_dz]
                                       # at the same density and radius (a 3-layer column)
    film_only: bool = True             # with film_dz: True = all algae in the film (same cells per m^2 as
                                       # uniform over dz_top), False = uniform over both sub-layers. Compare the
                                       # two at the same film_dz: splitting a layer is not exactly neutral in
                                       # the delta-Eddington adding-doubling solver (~0.005 albedo).

    def key(self):
        return (self.mode, self.rds_um, self.rho_top, self.rho_bottom, self.dz_top, self.dz_bottom, self.film_dz,
                self.film_only)

    @property
    def split(self) -> bool:
        return self.film_dz is not None and self.film_dz < self.dz_top - 1e-9


# Field counts are cells per mL of MELTWATER (melted surface ice, 1 mL = 1 g). BioSNICAR converts its
# cells/mL input to cells per kg of ice as conc / 917 * 1e6, i.e. per mL of SOLID ice
# (column_OPs.mix_in_impurities). Passing conc_field * 0.917 makes the column number of cells equal
# conc_field [cells/g] * rho * dz [g/m^2], as measured.
MELTWATER_TO_BIOSNICAR = 917.0 / 1000.0


def _layer_concs(spec: IceSpec, conc: float, unit: int):
    """Impurity concentration per model layer, with the meltwater-unit correction for cell counts and
    the areal-number-conserving rescaling for a thin algal film."""
    c = float(conc) * (MELTWATER_TO_BIOSNICAR if unit == 1 else 1.0)
    if spec.split:
        return [c * spec.dz_top / spec.film_dz, 0.0, 0.0] if spec.film_only else [c, c, 0.0]
    return [c, 0.0]


FIELD_SEASON_DOY = 196            # 15 July: centre of the S6 field campaigns used throughout


def sw_down_clear_sky(sza_deg: float, transmissivity: float | None = None, day_of_year: int = FIELD_SEASON_DOY) -> float:
    """Broadband clear-sky downwelling shortwave at the surface [W m^-2]:
    SW = S0 E0(doy) cos(SZA) T^(1/cos SZA) (Beer-Lambert bulk atmosphere, relative air mass
    1/cos SZA; E0 = 1 + 0.033 cos(2 pi doy / 365), the Earth-Sun distance factor used in the fit of T). Default T: fitted to PROMICE KAN_M clear-sky hours (0.919;
    empirical_data.clear_sky_transmissivity); pass measured fluxes via --sw-down instead."""
    if transmissivity is None:
        import empirical_data as ED
        transmissivity = ED.clear_sky_transmissivity()[0]
    mu = np.cos(np.radians(sza_deg))
    e0 = 1.0 + 0.033 * np.cos(2.0 * np.pi * day_of_year / 365.0)
    return float(S0 * e0 * mu * transmissivity ** (1.0 / mu))


class BioSNICARRunner:
    """Caches clean-ice optics per IceSpec and runs the solver for any impurity."""

    def __init__(self, biosnicar_root: str, incoming: int = 3, direct: int = 1,
                 band_nm=(300.0, 2500.0)):
        from biosnicar.drivers.setup_snicar import setup_snicar

        self.root = biosnicar_root
        self.input_file = os.path.join(biosnicar_root, "biosnicar", "inputs.yaml")
        (self._ice0, self._ill0, _rt, self.model_config, _pc, self._imps0) = setup_snicar(self.input_file)
        self.wvl_um = np.asarray(self.model_config.wavelengths, dtype=float)
        self.band = (self.wvl_um >= band_nm[0] / 1000.0) & (self.wvl_um <= band_nm[1] / 1000.0)
        self.incoming = incoming          # 3 = sub-Arctic summer (swnb_480bnd_sas)
        self.direct = direct
        self._ice_cache = {}
        self._ill_cache = {}

    # -- default BioSNICAR impurity (tier A) ------------------------------------
    def default_impurity(self, name: str = "glacier_algae"):
        imp = next(i for i in self._imps0 if i.name == name)
        return copy.deepcopy(imp)

    # -- ice column -------------------------------------------------------------
    def _base_ice(self, mode: str, n_layers: int = 2):
        """n-layer Ice object for a mode, refractive index computed once."""
        if ("base", mode, n_layers) not in self._ice_cache:
            ice = copy.deepcopy(self._ice0)
            lt = 0 if mode == "grains" else 1
            ice.nbr_lyr = n_layers
            ice.layer_type = [lt] * n_layers
            for attr in ("cdom", "shp", "water", "hex_side", "hex_length", "shp_fctr", "grain_ar", "lwc",
                         "rds", "rho", "dz"):
                v = getattr(ice, attr)
                setattr(ice, attr, (list(v) + [v[-1]] * n_layers)[:n_layers])
            ice.calculate_refractive_index(self.input_file)
            self._ice_cache[("base", mode, n_layers)] = ice
        return self._ice_cache[("base", mode, n_layers)]

    def available_radii(self, mode: str = "grains") -> np.ndarray:
        """Radii (um) tabulated in BioSNICAR's look-up table for this ice mode."""
        from biosnicar.optical_properties import column_OPs as cop

        if ("radii", mode) not in self._ice_cache:
            ice = self._base_ice(mode)
            path = cop._sphere_lut_path(self.model_config, ice) if mode == "grains" \
                else cop._bubbly_air_lut_path(self.model_config)
            self._ice_cache[("radii", mode)] = np.load(path)["radii"].astype(int)
        return self._ice_cache[("radii", mode)]

    def snap_radius(self, rds_um: float, mode: str = "grains") -> int:
        """Nearest tabulated radius (BioSNICAR's LUT steps are 20 um between 1 and 5 mm)."""
        r = self.available_radii(mode)
        return int(r[np.argmin(np.abs(r - rds_um))])

    def _fast_grain_ops(self, ice, rds: int):
        """Same numbers as get_layer_OPs for plain spherical grains (shp=0, no water
        coating), but indexing an in-memory copy of the LUT. BioSNICAR's LUT accessor
        decompresses the whole array from the .npz on every call (~40 ms each)."""
        from biosnicar.optical_properties import column_OPs as cop

        if any(s != 0 for s in ice.shp) or any(w > r for w, r in zip(ice.water, ice.rds)):
            return None
        key = ("lut", "grains")
        if key not in self._ice_cache:
            d = np.load(cop._sphere_lut_path(self.model_config, ice))
            self._ice_cache[key] = {k: d[k] for k in ("radii", "ss_alb", "ext_cff_mss", "asm_prm")}
            self._ice_cache[key]["idx"] = {int(r): i for i, r in enumerate(d["radii"])}
        lut = self._ice_cache[key]
        i = lut["idx"][int(rds)]
        n = ice.nbr_lyr
        return (np.tile(lut["ss_alb"][i], (n, 1)), np.tile(lut["asm_prm"][i], (n, 1)),
                np.tile(lut["ext_cff_mss"][i], (n, 1)))

    def ice(self, spec: IceSpec):
        """(ice object, layer SSA, g, mass extinction) for an IceSpec.

        Layer optics of granular ice depend only on grain radius (density enters later
        through the layer mass rho*dz), so they are cached per radius; bubbly ice also
        depends on density. The cheap Ice object copy carries this spec's rho/dz.
        """
        from biosnicar.optical_properties.column_OPs import get_layer_OPs

        if spec.key() in self._ice_cache:
            return self._ice_cache[spec.key()]
        rds = self.snap_radius(spec.rds_um, spec.mode)
        ice = copy.copy(self._base_ice(spec.mode))
        ice.rds = [rds, rds]
        ice.rho = [spec.rho_top, spec.rho_bottom]
        ice.dz = [spec.dz_top, spec.dz_bottom]
        op_key = ("ops", spec.mode, rds) + ((spec.rho_top, spec.rho_bottom) if spec.mode != "grains" else ())
        if op_key not in self._ice_cache:
            fast = self._fast_grain_ops(ice, rds) if spec.mode == "grains" else None
            self._ice_cache[op_key] = fast if fast is not None else get_layer_OPs(ice, self.model_config)
        ssa, g, mac = self._ice_cache[op_key]
        if spec.split:
            # same optics in the film and in the rest of the crust (same radius and density)
            ice = copy.copy(self._base_ice(spec.mode, 3))
            ice.rds = [rds, rds, rds]
            ice.rho = [spec.rho_top, spec.rho_top, spec.rho_bottom]
            ice.dz = [spec.film_dz, spec.dz_top - spec.film_dz, spec.dz_bottom]
            ssa, g, mac = (np.vstack([a[0], a[0], a[1]]) for a in (ssa, g, mac))
        out = (ice, ssa, g, mac)
        if len(self._ice_cache) < 20000:        # bounded memo for sweeps; MC samples rarely repeat
            self._ice_cache[spec.key()] = out
        return out

    @staticmethod
    def detach_luts():
        """Read BioSNICAR's cached look-up tables fully into memory. Its tables are lazy np.load
        archives sharing one file handle, which breaks when worker processes are forked."""
        from biosnicar.optical_properties import op_lookup
        for t in op_lookup._cache.values():
            if hasattr(t.data, "files"):
                t.data = {k: t.data[k] for k in t.data.files}

    def illumination(self, sza_deg: float):
        key = int(round(sza_deg))
        if key not in self._ill_cache:
            ill = copy.deepcopy(self._ill0)
            ill.solzen = key
            ill.incoming = self.incoming
            ill.direct = self.direct
            ill.calculate_irradiance()
            self._ill_cache[key] = ill
        return self._ill_cache[key]

    # -- one forward run ----------------------------------------------------------
    def run(self, spec: IceSpec, sza_deg: float, impurity=None, conc_cells_ml: float = 0.0):
        """Spectral albedo (480,) and normalised spectral irradiance (480,)."""
        from biosnicar.optical_properties.column_OPs import mix_in_impurities
        from biosnicar.rt_solvers.adding_doubling_solver import adding_doubling_solver

        ice, ssa_i, g_i, mac_i = self.ice(spec)
        imps = []
        if impurity is not None and conc_cells_ml > 0:
            impurity.conc = _layer_concs(spec, conc_cells_ml, impurity.unit)   # algae only at the surface
            imps = [impurity]
        tau, ssa, g, L = mix_in_impurities(ssa_i, g_i, mac_i, ice, imps, self.model_config)
        ill = self.illumination(sza_deg)
        out = adding_doubling_solver(tau, ssa, g, L, ice, ill, self.model_config)
        return np.asarray(out.albedo, dtype=float), np.asarray(ill.flx_slr, dtype=float), float(out.BBA)

    def run_multi(self, spec: IceSpec, sza_deg: float, impurities_concs):
        """Like run(), for several impurities at once: [(impurity, conc), ...] in the surface layer
        (e.g. two algal species, or algae + mineral dust). Units follow each impurity's .unit."""
        from biosnicar.optical_properties.column_OPs import mix_in_impurities
        from biosnicar.rt_solvers.adding_doubling_solver import adding_doubling_solver

        ice, ssa_i, g_i, mac_i = self.ice(spec)
        imps = []
        for imp, c in impurities_concs:
            if c > 0:
                imp.conc = _layer_concs(spec, c, imp.unit)
                imps.append(imp)
        tau, ssa, g, L = mix_in_impurities(ssa_i, g_i, mac_i, ice, imps, self.model_config)
        ill = self.illumination(sza_deg)
        out = adding_doubling_solver(tau, ssa, g, L, ice, ill, self.model_config)
        return np.asarray(out.albedo, dtype=float), np.asarray(ill.flx_slr, dtype=float), float(out.BBA)

    def broadband(self, albedo, flx):
        """Irradiance-weighted broadband albedo over the configured band (300-2500 nm)."""
        b = self.band
        return float(np.sum(albedo[b] * flx[b]) / np.sum(flx[b]))

    def forcing(self, albedo_clean, albedo, flx, sw_down: float):
        """Instantaneous surface shortwave forcing of the impurity [W m^-2] over the band:
        RF = SW_down * sum_lambda f(lambda) [alpha_clean - alpha], f = normalised spectrum."""
        b = self.band
        return float(sw_down * np.sum(flx[b] * (albedo_clean[b] - albedo[b])))


def export_lap_entry(path: str, stem: str, optics: dict):
    """Write one impurity in BioSNICAR's lap.npz key convention (for ADDING_DATA.md)."""
    np.savez_compressed(path, **{f"{stem}__ext_xsc": optics["ext_xsc"],
                                 f"{stem}__ss_alb": optics["ss_alb"],
                                 f"{stem}__asm_prm": optics["asm_prm"]})
