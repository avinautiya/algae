"""
Uncertain inputs of the molecule -> cell -> ice-surface chain, as probability
distributions, plus samplers.

Every parameter is a frozen scipy.stats distribution. Samplers work in the unit
hypercube u in (0,1)^D and map to physical values with the inverse CDF (ppf),
so Latin Hypercube (Monte Carlo) and Saltelli/Sobol' (sensitivity) designs use
exactly the same PDFs. Sobol' indices are invariant to this monotone transform.

Inputs are assumed independent (a requirement of the classical Sobol'
decomposition). If you later have data linking e.g. cell size and pigment
concentration, note that correlated inputs need Shapley effects instead.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "phase2"))

GROUPS = ("molecular", "cellular", "environmental")




@dataclass(frozen=True)
class Param:
    name: str
    label: str            # mathtext label for figures
    unit: str
    group: str            # molecular / cellular / environmental
    dist: object          # frozen scipy.stats distribution
    rationale: str
    transform: object = None  # applied after the inverse CDF (e.g. np.exp for a PDF defined in ln x)

    def ppf(self, u):
        v = self.dist.ppf(u)
        return v if self.transform is None else self.transform(v)


class _UnitDraw:
    """Uniform(0, 1) 'distribution' of a posterior-draw index (ppf is the identity)."""
    def ppf(self, u):
        return np.asarray(u, float)


def default_parameters(include_tier_d: bool = True, calibration: dict | None = None,
                       joint_calibration: bool = True) -> list[Param]:
    """Default PDFs, all from data. Molecular: marginal posteriors of the empirical TD-DFT calibration
    (phase2/tddft_calibration.py; pass `calibration` = Calibration.summary()). Cellular and
    environmental: published measurements (phase2/empirical_data.py, data/empirical/SOURCES.md).
    Molecular parameters: with joint_calibration (default) ONE input, cal_draw ~ U(0, 1), selects a row
    of the calibration's joint posterior samples (dE, FWHM, f, phi together), so their correlations
    (f-phi is strongly correlated) are kept and the Sobol' inputs stay independent; the molecular group
    index is then that of the joint draw. joint_calibration=False reproduces the former independent
    marginals (which dropped the correlations) for sensitivity."""
    if calibration is None:
        raise ValueError("default_parameters needs the TD-DFT calibration summary (tddft_calibration)")
    cal = calibration
    import empirical_data as ED
    _V = ED.s6_biovolume_um3()                                    # pooled, per-sample mean, SD, n
    _AR = [g["aspect"] for g in ED.species_geometry().values()]
    ci = ED.intracellular_concentration_kg_m3("phenolics")
    m, sd = ED.pigments_per_cell()["phenolics"]
    _CI = (ci, ci * sd / m)
    _B = ED.abundance_prior()
    _RHO = ED.ICE_DENSITY["weathering_crust"]
    _SSA = ED.ice_ssa_prior()
    _G = ED.phenolic_size_scaling()
    _TAU = ED.clear_sky_transmissivity()
    _DUST = ED.dust_prior()

    def _tn(lo, hi, m, sd):
        return stats.truncnorm((lo - m) / sd, (hi - m) / sd, loc=m, scale=sd)
    joint = [Param("cal_draw", r"TD-DFT calibration draw", "-", "molecular", _UnitDraw(),
                   "Index into the joint posterior samples of the TD-DFT calibration (dE, FWHM, f, phi), "
                   "phase2/tddft_calibration.py; correlations preserved.")]
    p = [
        # ---------------- molecular -------------------------------------------------
        Param("dE_ev", r"$\Delta E$ (TD-DFT shift)", "eV", "molecular",
              _tn(-1.0, 1.0, cal["dE"]["mean"], cal["dE"]["sd"]),
              f"Band shift of the Phase 1 spectrum, N({cal['dE']['mean']:+.3f}, {cal['dE']['sd']:.3f}) eV: posterior of "
              "the fit to HPLC spectra of the isolated pigment (Williamson et al. 2020)."),
        Param("f_scale", r"$f$ scale factor", "-", "molecular",
              stats.lognorm(s=cal["log_f"]["sd"], scale=np.exp(cal["log_f"]["mean"])),
              f"Oscillator-strength scale incl. glucoside -> phenol-equivalent units, log-normal (median "
              f"{np.exp(cal['log_f']['mean']):.3g}, sigma_ln {cal['log_f']['sd']:.3f}): posterior of the fit to the "
              "measured S6 extract MAC (Williamson et al. 2020)."),
        Param("fwhm_ev", r"Band FWHM", "eV", "molecular",
              _tn(0.05, 2.0, cal["w"]["mean"], cal["w"]["sd"]),
              f"Gaussian band FWHM, N({cal['w']['mean']:.3f}, {cal['w']['sd']:.3f}) eV: posterior of the HPLC-shape fit."),
        # ---------------- cellular (measured; phase2/empirical_data.py, data/empirical/SOURCES.md) ----
        Param("cell_volume_um3", r"Cell volume $V$", r"$\mu$m$^3$", "cellular",
              stats.truncnorm((400.0 - _V[1]) / _V[2], np.inf, loc=_V[1], scale=_V[2]),
              f"Mean glacier-algal biovolume per cell, N({_V[1]:.0f}, {_V[2]:.0f}) um^3: mean and SD over "
              f"{_V[3]} S6 surface-ice samples (Williamson et al. 2020 counts), truncated at 400 um^3."),
        Param("cell_aspect", r"Aspect $L/d$", "-", "cellular",
              stats.uniform(min(_AR), max(_AR) - min(_AR)),
              f"Cylinder length/diameter, U({min(_AR):.2f}, {max(_AR):.2f}): between the mean ratios of "
              "A. alaskanum and A. nordenskioeldii populations (Prochazkova et al. 2021, Table 2)."),
        Param("c_internal", r"Pigment conc. $c_i$", r"kg m$^{-3}$", "cellular",
              stats.truncnorm((0.5 - _CI[0]) / _CI[1], np.inf, loc=_CI[0], scale=_CI[1]),
              f"Intracellular phenolic concentration N({_CI[0]:.1f}, {_CI[1]:.1f}) kg m^-3: phenolics per cell "
              "(mean, SD over 53 samples) / pooled S6 biovolume (Williamson et al. 2020); truncated at 0.5."),
        Param("conc_size_exponent", r"Conc. size exponent $\gamma$", "-", "cellular",
              stats.norm(_G["gamma"], _G["gamma_se"]),
              f"Intracellular concentration ~ V^gamma, gamma ~ N({_G['gamma']:.2f}, {_G['gamma_se']:.2f}) (jackknife SE): "
              f"maximum-likelihood fit of phenolics per cell vs biovolume per cell over {_G['n']} S6 samples "
              "(Williamson et al. 2020); gamma = 0 is equal concentration in both species."),
        # ---------------- environmental (measured) ----------------------------------
        Param("ice_ssa", r"Ice SSA", r"m$^2$ kg$^{-1}$", "environmental",
              stats.lognorm(s=_SSA[1], scale=np.exp(_SSA[0])),
              f"Specific surface area of bubbly bare ice, ln SSA ~ N({_SSA[0]:.3f}, {_SSA[1]:.3f}): {_SSA[2]} measurements "
              "(Cooper et al. 2021 Greenland; Dadic et al. 2013 micro-CT), studies weighted equally; converted to "
              "BioSNICAR's bubble radius at the bottom-layer density."),
        Param("rho_top", r"Surface density $\rho$", r"kg m$^{-3}$", "environmental",
              stats.uniform(_RHO["lo"], _RHO["hi"] - _RHO["lo"]),
              f"Weathering-crust density, U({_RHO['lo']:.0f}, {_RHO['hi']:.0f}) kg m^-3: measured range at S6 "
              "(Cooper et al. 2018)."),
        Param("conc_cells_ml", r"Cell abundance", r"cells mL$^{-1}$", "environmental",
              stats.lognorm(s=_B[1] * np.log(10.0), scale=10.0 ** _B[0]),
              f"Algal abundance, log10 B ~ N({_B[0]:.2f}, {_B[1]:.2f}): {_B[2]} S6 surface-ice samples with "
              "cells > 0 (Williamson et al. 2020 counts)."),
        Param("dust_ppb", r"Mineral dust", "ppb", "environmental",
              stats.lognorm(s=_DUST[1], scale=np.exp(_DUST[0])),
              f"Surface mineral-dust mass mixing ratio, log-normal (median {np.exp(_DUST[0]) / 1e3:.0f} ug/g, "
              f"sigma_ln {_DUST[1]:.2f}), moment-matched to the measured S6 loading 342 ug/g, relative SD 0.49 "
              "(Cook et al. 2020)."),
        Param("transmissivity", r"Clear-sky $T$", "-", "environmental",
              _tn(0.5, 1.0, _TAU[0], _TAU[1]),
              f"Bulk clear-sky transmissivity N({_TAU[0]:.3f}, {_TAU[1]:.3f}): {_TAU[2]} clear-sky hours at PROMICE KAN_M "
              "(June-August, all years)."),
    ]
    if joint_calibration:
        p = joint + [q for q in p if q.name not in ("dE_ev", "f_scale", "fwhm_ev")]
        return p
    if include_tier_d:
        p += [
            Param("fe_fraction", r"Fe-complexed fraction $\phi$", "-", "molecular",
                  _tn(0.0, 1.0, cal["phi"]["mean"], max(cal["phi"]["sd"], 1e-3)),
                  f"Fraction of pigment complexed with Fe, N({cal['phi']['mean']:.3f}, {cal['phi']['sd']:.3f}) on [0, 1]: "
                  "posterior of the fit of the measured Fe-purpurogallin increment (Prochazkova et al. 2025) to the "
                  "S6 extract MAC (tier D only)."),
        ]
    return p


class ParameterSpace:
    def __init__(self, params: list[Param]):
        self.params = list(params)
        self.names = [p.name for p in self.params]
        self.D = len(self.params)

    # ---- transforms -------------------------------------------------------------
    def from_unit(self, U: np.ndarray) -> np.ndarray:
        U = np.clip(np.asarray(U, dtype=float), 1e-9, 1.0 - 1e-9)
        return np.column_stack([p.ppf(U[:, j]) for j, p in enumerate(self.params)])

    def as_dicts(self, X: np.ndarray) -> list[dict]:
        return [dict(zip(self.names, row)) for row in X]

    def medians(self) -> dict:
        return {p.name: float(p.ppf(0.5)) for p in self.params}

    def quantiles(self, q: float) -> dict:
        return {p.name: float(p.ppf(q)) for p in self.params}

    # ---- designs ----------------------------------------------------------------
    def lhs(self, n: int, seed: int = 2024) -> np.ndarray:
        """Latin Hypercube sample (n x D) in physical units (scipy.stats.qmc, optimised
        with random-CD to reduce spurious correlations between columns)."""
        from scipy.stats import qmc

        sampler = qmc.LatinHypercube(d=self.D, seed=seed, optimization="random-cd")
        return self.from_unit(sampler.random(n))

    def salib_problem(self, groups: bool = False) -> dict:
        prob = {"num_vars": self.D, "names": self.names, "bounds": [[0.0, 1.0]] * self.D}
        if groups:
            prob["groups"] = [p.group for p in self.params]
        return prob

    def saltelli(self, n_base: int, second_order: bool = True, groups: bool = False, seed: int = 2024):
        """Saltelli (2010) design via SALib (scrambled Sobol' sequence), in physical units.
        n_base must be a power of two. Returns (X, problem)."""
        from SALib.sample import sobol as sobol_sample

        if n_base & (n_base - 1):
            raise ValueError("n_base must be a power of two for Sobol' sequences")
        prob = self.salib_problem(groups)
        U = sobol_sample.sample(prob, n_base, calc_second_order=second_order, seed=seed)
        return self.from_unit(U), prob

    def table(self):
        import pandas as pd
        rows = []
        for p in self.params:
            rows.append(dict(parameter=p.name, group=p.group, unit=p.unit, median=float(p.ppf(0.5)),
                             p2_5=float(p.ppf(0.025)), p97_5=float(p.ppf(0.975)), rationale=p.rationale))
        return pd.DataFrame(rows)
