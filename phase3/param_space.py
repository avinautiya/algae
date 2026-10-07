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

GROUPS = ("molecular", "cellular", "environmental")


@dataclass(frozen=True)
class Param:
    name: str
    label: str            # mathtext label for figures
    unit: str
    group: str            # molecular / cellular / environmental
    dist: object          # frozen scipy.stats distribution
    rationale: str

    def ppf(self, u):
        return self.dist.ppf(u)


def default_parameters(include_tier_d: bool = True) -> list[Param]:
    """Default PDFs. Edit here (or pass your own list) to change the uncertainty model."""
    p = [
        # ---------------- molecular -------------------------------------------------
        Param("dE_ev", r"$\Delta E$ (TD-DFT shift)", "eV", "molecular",
              stats.truncnorm(-2.0, 2.0, loc=0.0, scale=0.075),
              "Systematic TD-DFT excitation-energy error: N(0, 0.075 eV) truncated at +/-0.15 eV "
              "(typical B3LYP/CAM-B3LYP error for pi->pi* bands of polyphenols)."),
        Param("f_scale", r"$f$ scale factor", "-", "molecular",
              stats.lognorm(s=0.20, scale=1.0),
              "Multiplicative oscillator-strength error, log-normal with median 1 and "
              "sigma_ln = 0.20 (~+/-20 %, 1 sigma)."),
        Param("fwhm_ev", r"Band FWHM", "eV", "molecular",
              stats.uniform(0.25, 0.15),
              "Vibronic + inhomogeneous solvent broadening, U(0.25, 0.40) eV around the "
              "Phase 1 choice of 0.30 eV."),
        # ---------------- cellular --------------------------------------------------
        Param("cell_length_um", r"Cell length $L$", r"$\mu$m", "cellular",
              stats.uniform(10.0, 20.0), "Ancylonema cell length, U(10, 30) um."),
        Param("cell_diameter_um", r"Cell diameter $d$", r"$\mu$m", "cellular",
              stats.uniform(5.0, 10.0), "Ancylonema cell diameter, U(5, 15) um."),
        Param("c_internal", r"Pigment conc. $c_i$", r"kg m$^{-3}$", "cellular",
              stats.loguniform(10.0, 200.0),
              "Intracellular phenolic concentration, log-uniform 10-200 kg m^-3 "
              "(order-of-magnitude prior; replace with HPLC-based per-cell estimates)."),
        # ---------------- environmental ---------------------------------------------
        Param("grain_um", r"Ice grain radius", r"$\mu$m", "environmental",
              stats.uniform(1000.0, 2000.0), "Weathering-crust grain radius, U(1, 3) mm."),
        Param("rho_top", r"Surface density $\rho$", r"kg m$^{-3}$", "environmental",
              stats.uniform(500.0, 300.0), "Weathering-crust density, U(500, 800) kg m^-3."),
        Param("conc_cells_ml", r"Cell abundance", r"cells mL$^{-1}$", "environmental",
              stats.loguniform(1e3, 1e5), "Algal abundance, log-uniform 10^3-10^5 cells mL^-1."),
    ]
    if include_tier_d:
        p += [
            Param("lmct_eps", r"LMCT $\varepsilon_{max}$", r"M$^{-1}$cm$^{-1}$", "molecular",
                  stats.uniform(3000.0, 2000.0),
                  "Fe(III)<-phenolate LMCT molar absorptivity per Fe, U(3000, 5000) (tier D only)."),
            Param("lmct_center_nm", r"LMCT $\lambda_{max}$", "nm", "molecular",
                  stats.uniform(520.0, 100.0),
                  "LMCT band centre, U(520, 620) nm (tier D only)."),
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
        return {p.name: float(p.dist.median()) for p in self.params}

    def quantiles(self, q: float) -> dict:
        return {p.name: float(p.dist.ppf(q)) for p in self.params}

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
            d = p.dist
            rows.append(dict(parameter=p.name, group=p.group, unit=p.unit, median=d.median(),
                             p2_5=d.ppf(0.025), p97_5=d.ppf(0.975), rationale=p.rationale))
        return pd.DataFrame(rows)
