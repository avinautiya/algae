"""
Fixed-optics forward model for the baseline inversion: BioSNICAR (via phase2/biosnicar_bridge) with the
shipped empirical optics of glacier algae, Greenland dust and black carbon, on the bubbly-ice column of
phase4/emulator.EmulatorConfig.

The bridge converts meltwater cell counts to BioSNICAR units internally (x 0.917), so this module passes RAW
meltwater counts. The optical radius is continuous: albedo is interpolated linearly in ln(radius) between the two
tabulated bubbly-ice radii that bracket it, so the objective is smooth for a gradient-based fit.
"""

from __future__ import annotations

import copy
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))
sys.path.insert(0, HERE)


class FixedOpticsForward:
    def __init__(self, biosnicar_root: str, cfg: dict):
        import biosnicar_bridge as BB
        self.BB = BB
        self.root = BB.locate_biosnicar(biosnicar_root)
        self.runner = BB.BioSNICARRunner(self.root)
        ic = cfg["forward"]["ice_column"]
        self.ice = dict(rho_top=ic["rho_crust_kg_m3"], rho_bottom=ic["rho_ice_kg_m3"], dz_top=ic["dz_crust_m"],
                        dz_bottom=ic["dz_ice_m"])
        self.algae = self.runner.default_impurity("glacier_algae")
        self.dust = self.runner.default_impurity("dust")
        self.bc = self.runner.default_impurity("black_carbon")
        self.wl_nm = self.runner.wvl_um * 1000.0
        self.radii = np.asarray(self.runner.available_radii("bubbly"), dtype=float)
        self._flux = {}
        assert (self.algae.unit, self.dust.unit, self.bc.unit) == (1, 0, 0), "unit conventions changed in BioSNICAR"

    def _run(self, r_um: float, sza: float, cells: float, dust_ppb: float, bc_ppb: float) -> np.ndarray:
        spec = self.BB.IceSpec(rds_um=float(r_um), mode="bubbly", **self.ice)
        a, _flx, _bba = self.runner.run_multi(spec, sza, [(self.algae, cells), (self.dust, dust_ppb), (self.bc, bc_ppb)])
        return a

    def albedo(self, log10_cells: float, log10_dust_ppb: float, log10_bc_ppb: float, ln_r_um: float, sza: float):
        """Spectral albedo on the 480-band model grid (hemispherical, direct beam at `sza`)."""
        cells, dust, bc = 10.0 ** log10_cells, 10.0 ** log10_dust_ppb, 10.0 ** log10_bc_ppb
        r = float(np.clip(np.exp(ln_r_um), self.radii[0], self.radii[-1]))
        j = int(np.clip(np.searchsorted(self.radii, r), 1, len(self.radii) - 1))
        r0, r1 = self.radii[j - 1], self.radii[j]
        w = (np.log(r) - np.log(r0)) / (np.log(r1) - np.log(r0))
        a0 = self._run(r0, sza, cells, dust, bc)
        if w <= 1e-9:
            return a0
        a1 = self._run(r1, sza, cells, dust, bc)
        return (1.0 - w) * a0 + w * a1

    def flux(self, sza: float) -> np.ndarray:
        key = int(round(sza))
        if key not in self._flux:
            self._flux[key] = np.asarray(self.runner.illumination(sza).flx_slr, dtype=float)
        return self._flux[key]
