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

# ln(bubbly-ice radius / um): mean and SD of the posterior-mean radii retrieved from the 31 Cook et al.
# (2020) field spectra by phase4/field_validation.py (empirical configuration, Williamson MACs).
FIELD_LN_RADIUS = (7.386, 0.977)


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


def default_parameters(include_tier_d: bool = True) -> list[Param]:
    """Default PDFs. Cellular and environmental PDFs come from published measurements
    (phase2/empirical_data.py); the molecular PDFs are user-specified TD-DFT error models."""
    import empirical_data as ED
    _V = ED.s6_biovolume_um3()                                    # pooled, per-sample mean, SD, n
    _AR = [g["aspect"] for g in ED.species_geometry().values()]
    ci = ED.intracellular_concentration_kg_m3("phenolics")
    m, sd = ED.pigments_per_cell()["phenolics"]
    _CI = (ci, ci * sd / m)
    _B = ED.abundance_prior()
    _RHO = ED.ICE_DENSITY["weathering_crust"]
    _R_LO, _R_HI = ED.ICE_RADIUS_BOUNDS_UM
    _LR = FIELD_LN_RADIUS
    p = [
        # ---------------- molecular -------------------------------------------------
        Param("dE_ev", r"$\Delta E$ (TD-DFT shift)", "eV", "molecular",
              stats.truncnorm(-2.0, 2.0, loc=0.0, scale=0.075),
              "Systematic TD-DFT excitation-energy error: N(0, 0.075 eV) truncated at +/-0.15 eV "
              "(USER-SPECIFIED, not measured: no experimental spectrum of the glucoside to calibrate against)."),
        Param("f_scale", r"$f$ scale factor", "-", "molecular",
              stats.lognorm(s=0.20, scale=1.0),
              "Multiplicative oscillator-strength error, log-normal with median 1 and "
              "sigma_ln = 0.20 (~+/-20 %, 1 sigma) (USER-SPECIFIED, not measured)."),
        Param("fwhm_ev", r"Band FWHM", "eV", "molecular",
              stats.uniform(0.25, 0.15),
              "Vibronic + inhomogeneous solvent broadening, U(0.25, 0.40) eV around the "
              "Phase 1 choice of 0.30 eV (USER-SPECIFIED, not measured)."),
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
        # ---------------- environmental (measured) ----------------------------------
        Param("grain_um", r"Bubbly-ice radius", r"$\mu$m", "environmental",
              stats.truncnorm((np.log(_R_LO) - _LR[0]) / _LR[1], (np.log(_R_HI) - _LR[0]) / _LR[1],
                              loc=_LR[0], scale=_LR[1]),
              f"Bubbly-ice optical radius, log-normal (ln r ~ N({_LR[0]:.2f}, {_LR[1]:.2f})) on "
              f"[{_R_LO:.0f}, {_R_HI:.0f}] um: fitted to the radii retrieved by the Phase 4 inversion from the "
              "31 Cook et al. (2020) field spectra.", transform=np.exp),
        Param("rho_top", r"Surface density $\rho$", r"kg m$^{-3}$", "environmental",
              stats.uniform(_RHO["lo"], _RHO["hi"] - _RHO["lo"]),
              f"Weathering-crust density, U({_RHO['lo']:.0f}, {_RHO['hi']:.0f}) kg m^-3: measured range at S6 "
              "(Cooper et al. 2018)."),
        Param("conc_cells_ml", r"Cell abundance", r"cells mL$^{-1}$", "environmental",
              stats.lognorm(s=_B[1] * np.log(10.0), scale=10.0 ** _B[0]),
              f"Algal abundance, log10 B ~ N({_B[0]:.2f}, {_B[1]:.2f}): {_B[2]} S6 surface-ice samples with "
              "cells > 0 (Williamson et al. 2020 counts)."),
    ]
    if include_tier_d:
        p += [
            Param("lmct_eps", r"LMCT $\varepsilon_{max}$", r"M$^{-1}$cm$^{-1}$", "molecular",
                  stats.uniform(3000.0, 2000.0),
                  "Fe(III)<-phenolate LMCT molar absorptivity per Fe, U(3000, 5000) (tier D only; PROVISIONAL "
                  "surrogate - no measured spectrum of the algal Fe-phenolic complex is published)."),
            Param("lmct_center_nm", r"LMCT $\lambda_{max}$", "nm", "molecular",
                  stats.uniform(520.0, 100.0),
                  "LMCT band centre, U(520, 620) nm (tier D only; PROVISIONAL, as above)."),
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
