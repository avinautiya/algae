"""Additive adapter for existing Phase 4 cell optics; no production-pipeline edits.

Uses the private Phase 4 builder deliberately; its shape/API is regression-tested
without claiming that mock tests validate BioSNICAR or pigment calibration.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

from .coupling import spectral_budget, _array


class AlgaeOpticsAdapter:
    """Paired no-algae/algae spectra with identical dust, ice and illumination.

    Construct a separate adapter per process and illumination configuration.
    direct=1/0 selects direct/diffuse in the installed BioSNICAR implementation.
    This is offline coupling, not an E3SM/MAR source-code integration.
    """

    def __init__(self, cfg, *, phase1_l2=None, biosnicar=None, demo=False,
                 direct=1, builder_factory=None):
        if direct not in (0, 1):
            raise ValueError("direct must be 0 (diffuse) or 1 (direct)")
        if builder_factory is None:
            sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "phase4"))
            from emulator import _Builder
            builder_factory = _Builder
        self.builder = builder_factory(cfg, phase1_l2, demo, biosnicar)
        self.builder.runner.direct = direct
        self.builder.runner._ill_cache.clear()

    def evaluate(self, *, cells_ml_meltwater, f_n, r_um, dust_ppb,
                 sza_deg, sw_down_w_m2):
        """One state/illumination; output bin flux and paired albedos for a host.

        Supply measured/host broadband SW; never normalize modeled RF using a
        measured overpass flux. Spectral distribution follows BioSNICAR incoming.
        Cloud-induced spectral changes need appropriate host spectral irradiance,
        which can instead be passed to coupling.spectral_budget directly.
        """
        B = float(_array(cells_ml_meltwater, "cells_ml_meltwater", 0))
        fn = float(_array(f_n, "f_n", 0, 1))
        radius = float(_array(r_um, "r_um", 0))
        dust = float(_array(dust_ppb, "dust_ppb", 0))
        zenith = float(_array(sza_deg, "sza_deg", 0, 180))
        sw = float(_array(sw_down_w_m2, "sw_down_w_m2", 0))
        if radius <= 0:
            raise ValueError("r_um must be positive")
        if zenith >= 90:
            if sw != 0:
                raise ValueError("nighttime nonzero SW requires a separately specified diffuse treatment")
            raise ValueError("skip nighttime RT and use zero absorbed SW for that interval")
        b = self.builder
        spec = b.spec(radius)
        # Reconciled with claude/sweet-noether-6pryo4: biosnicar_bridge._layer_concs places cell counts
        # (unit 1) in the film when film_only, and keeps mass-based dust (unit 0) uniform over the crust,
        # so film-only algae plus dust is a defined configuration (dust profile = uniform, stated).
        background = [(b.dust, dust)] if dust > 0 else []
        if b.cfg.model == "ours":
            algae = [(b.imps["nordenskioeldii"], B * fn),
                     (b.imps["alaskanum"], B * (1 - fn))]
        elif b.cfg.model == "tierA":
            algae = [(b.imps["tierA"], B)]
        else:
            raise ValueError("adapter supports ours or tierA")
        a0, weights, _ = b.runner.run_multi(spec, zenith, background)
        a1, weights1, _ = b.runner.run_multi(spec, zenith, background + algae)
        weights = _array(weights, "spectral illumination", 0)
        if weights.sum() <= 0 or not np.allclose(weights, weights1, rtol=1e-10, atol=1e-12):
            raise ValueError("paired illumination must match and have positive total weight")
        # Match the same 300–2500 nm band used in existing broadband diagnostics.
        mask = b.runner.band
        if weights[mask].sum() <= 0:
            raise ValueError("selected spectral band has no illumination")
        flux = weights[mask] / weights[mask].sum() * sw
        budget = spectral_budget(np.asarray(a0)[mask], np.asarray(a1)[mask], flux)
        return dict(wavelength_nm=np.asarray(b.runner.wvl_um)[mask] * 1000,
                    albedo_without_algae=np.asarray(a0)[mask],
                    albedo_with_algae=np.asarray(a1)[mask], flux_bin_w_m2=flux,
                    spectral_band_nm=[300, 2500],
                    broadband_sw_interpretation="all supplied SW allocated to 300–2500 nm; band-limited approximation",
                    **budget)
