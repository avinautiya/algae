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
    """Two-layer bare-ice column: algae-bearing weathering crust over clean ice.

    mode='grains': granular ice spheres of radius `rds_um` (layer_type 0)
    mode='bubbly': solid ice with Fresnel surface, `rds_um` = bubble radius (layer_type 1)
    """
    rds_um: float = 1500.0
    rho_top: float = 650.0
    rho_bottom: float = 850.0
    dz_top: float = 0.02
    dz_bottom: float = 2.0
    mode: str = "grains"

    def key(self):
        return (self.mode, self.rds_um, self.rho_top, self.rho_bottom, self.dz_top, self.dz_bottom)


def sw_down_clear_sky(sza_deg: float, transmissivity: float = 0.75) -> float:
    """Broadband clear-sky downwelling shortwave at the surface [W m^-2]:
    SW = S0 cos(SZA) T^(1/cos SZA) (Beer-Lambert bulk atmosphere, relative air mass
    1/cos SZA). Default T = 0.75 gives ~560 W m^-2 at SZA 50 deg, consistent with
    clear-sky ablation-zone observations; pass measured fluxes (e.g. PROMICE AWS)
    via --sw-down to replace it."""
    mu = np.cos(np.radians(sza_deg))
    return float(S0 * mu * transmissivity ** (1.0 / mu))


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
    def ice(self, spec: IceSpec):
        from biosnicar.optical_properties.column_OPs import get_layer_OPs

        if spec.key() not in self._ice_cache:
            ice = copy.deepcopy(self._ice0)
            lt = 0 if spec.mode == "grains" else 1
            ice.nbr_lyr = 2
            ice.layer_type = [lt, lt]
            ice.rds = [int(spec.rds_um)] * 2
            ice.rho = [spec.rho_top, spec.rho_bottom]
            ice.dz = [spec.dz_top, spec.dz_bottom]
            for attr in ("cdom", "shp", "water", "hex_side", "hex_length", "shp_fctr", "grain_ar", "lwc"):
                v = getattr(ice, attr)
                setattr(ice, attr, (list(v) + [v[-1]] * 2)[:2])
            ice.calculate_refractive_index(self.input_file)
            ssa, g, mac = get_layer_OPs(ice, self.model_config)
            self._ice_cache[spec.key()] = (ice, ssa, g, mac)
        return self._ice_cache[spec.key()]

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
            impurity.conc = [float(conc_cells_ml), 0.0]       # algae only in the surface layer
            imps = [impurity]
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
