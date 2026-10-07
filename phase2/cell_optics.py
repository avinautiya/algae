"""
Molecular MAC (Phase 1)  ->  single-cell optical properties on BioSNICAR's 480-band grid.

BioSNICAR treats algae given in cells/mL with three spectral arrays per impurity
(classes/impurity.py, column_OPs.mix_in_impurities):

    ext_xsc  extinction cross-section per cell         [m^2 cell^-1]
    ss_alb   single-scattering albedo                  [-]
    asm_prm  asymmetry parameter g                     [-]

This module builds those arrays for the model tiers from first principles:

  * absorption: pigment MAC x pigment mass per cell, either unpackaged (tier B)
    or with the Duysens packaging factor Q* (tiers C and D; packaging.py),
    plus the (packaged) absorption of intracellular water;
  * extinction: geometric-optics limit, sigma_ext = 2 <A_proj> (extinction
    paradox; the cells have size parameters 2 pi r / lambda ~ 40-600);
  * scattering: sigma_sca = sigma_ext - sigma_abs(packaged). Tier B keeps the
    SAME scattering as tier C and only replaces the absorption, so B vs C
    isolates the packaging effect;
  * asymmetry parameter: g = 0.96 for every tier by default (BioSNICAR's
    empirical glacier-algae value, so tiers differ only in absorption), or the
    van Diedenhoven et al. (2014) geometric-optics parameterisation BioSNICAR
    uses for algal cylinders (g_mode='vd2014'). The vd2014 single-scattering
    albedo is always returned as an independent cross-check of our
    packaging-based SSA.

Wavelength handling (task spec): the molecular MAC is used inside a window
(default 350-800 nm). Above the window pigment absorption is set to zero
(phenolics do not absorb in the NIR). Below it the MAC is held at its value at
the window edge ('hold'), or optionally taken from the TD-DFT spectrum down to
its lowest computed wavelength ('molecular') or set to zero ('zero').
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd

from pigment_packaging import CellGeometry, q_star

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase1"))
import spectra as p1spectra  # noqa: E402  (Phase 1 broadening utilities, reused verbatim)

HC_EV_NM = p1spectra.HC_EV_NM
WVL_480_UM = np.round(np.arange(0.205, 4.999, 0.01), 3)      # BioSNICAR band centres
WVL_480_NM = WVL_480_UM * 1000.0


# --------------------------------------------------------------------------- #
# Phase 1 molecular spectra                                                    #
# --------------------------------------------------------------------------- #
@dataclass
class MolecularSpectrum:
    """TD-DFT stick spectrum (preferred) or a continuous MAC curve."""
    name: str
    molar_mass: float                     # g/mol
    fwhm_ev: float = 0.3
    energies_ev: np.ndarray | None = None
    osc: np.ndarray | None = None
    wl_nm: np.ndarray | None = None       # continuous fallback
    mac: np.ndarray | None = None         # Napierian, m^2 kg^-1
    source: str = ""

    @property
    def has_sticks(self) -> bool:
        return self.energies_ev is not None

    def mac_at(self, wl_nm, shift_ev: float = 0.0, extra_fwhm_ev: float = 0.0):
        """Napierian MAC [m^2 kg^-1] at `wl_nm`.

        shift_ev < 0 red-shifts every band; extra_fwhm_ev adds Gaussian broadening
        in quadrature (used for the aggregated tier D).
        """
        wl = np.asarray(wl_nm, dtype=float)
        e = HC_EV_NM / wl
        fwhm = float(np.hypot(self.fwhm_ev, extra_fwhm_ev))
        if self.has_sticks:
            eps = p1spectra.gaussian_broaden(e, self.energies_ev + shift_ev, self.osc, fwhm)
            return p1spectra.epsilon_to_mac(eps, self.molar_mass, napierian=True)
        # continuous curve: shift/broaden numerically on a uniform energy grid
        e_src = HC_EV_NM / self.wl_nm[::-1]
        m_src = self.mac[::-1]
        if extra_fwhm_ev > 0:
            de = 0.002
            eg = np.arange(e_src.min(), e_src.max(), de)
            mg = np.interp(eg, e_src, m_src)
            sig = extra_fwhm_ev * p1spectra.FWHM_TO_SIGMA
            k = np.arange(-int(4 * sig / de), int(4 * sig / de) + 1) * de
            ker = np.exp(-0.5 * (k / sig) ** 2)
            mg = np.convolve(mg, ker / ker.sum(), mode="same")
            e_src, m_src = eg, mg
        return np.interp(e - shift_ev, e_src, m_src, left=0.0, right=0.0)


def load_phase1(phase1_dir: str, level: str, functional: str = "B3LYP") -> MolecularSpectrum:
    """Read a Phase 1 output folder (results/<level>/ from phase1/run_phase1.py)."""
    stem = os.path.join(phase1_dir, f"{level}_{functional}")
    with open(os.path.join(phase1_dir, "summary.json")) as fh:
        summ = json.load(fh)
    states = pd.read_csv(stem + "_states.csv")
    napierian = summ.get("napierian", True)
    if not napierian:
        raise ValueError("Phase 1 run used --decadic; rerun Phase 1 with the default Napierian MAC")
    spec = pd.read_csv(stem + "_spectrum.csv")
    return MolecularSpectrum(
        name=summ.get("name", level), molar_mass=float(summ["molar_mass_g_mol"]),
        fwhm_ev=float(summ.get("fwhm_ev", 0.3)),
        energies_ev=states["Energy_eV"].to_numpy(), osc=states["Oscillator_Strength"].to_numpy(),
        wl_nm=spec["Wavelength_nm"].to_numpy(), mac=spec["MAC_estimated"].to_numpy(),
        source=f"Phase 1 TD-{functional} ({summ.get('basis', '?')}, {summ.get('solvent', '?')})",
    )


def load_mac_csv(path: str, name: str, molar_mass: float) -> MolecularSpectrum:
    """Any continuous MAC table with columns Wavelength_nm, MAC_estimated (Napierian m^2/kg).
    Use this to plug in the Level 3/4 (Fe-complex / aggregate) spectra once computed."""
    df = pd.read_csv(path)
    return MolecularSpectrum(name=name, molar_mass=molar_mass, wl_nm=df["Wavelength_nm"].to_numpy(),
                             mac=df["MAC_estimated"].to_numpy(), source=os.path.basename(path))


def demo_spectrum(biosnicar_root: str) -> MolecularSpectrum:
    """DEMO ONLY: BioSNICAR's tabulated purpurogallin-type phenol MAC (data/pigments/ppg.csv,
    m^2 mg^-1 on a 1-nm grid from 200 nm). Lets the pipeline run before Phase 1 results
    exist; every output produced from it is labelled DEMO."""
    v = pd.read_csv(os.path.join(biosnicar_root, "data", "pigments", "ppg.csv"),
                    header=None, encoding="utf-8-sig").iloc[:, -1].to_numpy(dtype=float)
    wl = 200.0 + np.arange(v.size)
    keep = (wl >= 250) & (wl <= 900)
    return MolecularSpectrum(name="DEMO ppg (BioSNICAR)", molar_mass=426.33, wl_nm=wl[keep],
                             mac=v[keep] * 1e6, source="DEMO: BioSNICAR data/pigments/ppg.csv")


# --------------------------------------------------------------------------- #
# Tier D: Fe(III)-phenolic complexation + aggregation                           #
# --------------------------------------------------------------------------- #
@dataclass
class FePhenolicSurrogate:
    """PROVISIONAL stand-in for the Level 3/4 (Fe-complexed, aggregated) MAC.

    Replace with TD-DFT Level 3/4 output via `load_mac_csv` as soon as it exists.
    Parameters follow typical literature ranges for Fe(III)-catecholate /
    Fe(III)-galloyl chromophores (ligand-to-metal charge transfer, LMCT):
      * LMCT band centred ~500-600 nm (bis-catecholate ~570 nm), molar absorptivity
        ~3000-5000 M^-1 cm^-1 per Fe (decadic);
      * aggregation (pi-stacking) modelled as a small red shift and extra
        inhomogeneous broadening of the ligand pi->pi* bands.
    All values are configurable and must be justified/varied in the paper.
    """
    lmct_center_nm: float = 570.0
    lmct_fwhm_ev: float = 0.70
    lmct_eps_per_fe: float = 4000.0      # L mol^-1 cm^-1, decadic, per Fe
    fe_per_ligand: float = 0.5           # bis-complex: 1 Fe per 2 ligands
    agg_shift_ev: float = -0.10          # negative = red shift
    agg_extra_fwhm_ev: float = 0.25

    def mac(self, ligand: MolecularSpectrum, wl_nm):
        """MAC per kg of LIGAND (same molecule count per cell as tiers B/C)."""
        wl = np.asarray(wl_nm, dtype=float)
        base = ligand.mac_at(wl, shift_ev=self.agg_shift_ev, extra_fwhm_ev=self.agg_extra_fwhm_ev)
        e = HC_EV_NM / wl
        e0 = HC_EV_NM / self.lmct_center_nm
        sig = self.lmct_fwhm_ev * p1spectra.FWHM_TO_SIGMA
        eps = self.lmct_eps_per_fe * self.fe_per_ligand * np.exp(-0.5 * ((e - e0) / sig) ** 2)
        lmct = np.log(10.0) * 0.1 * eps / ligand.molar_mass * 1000.0     # m^2 kg^-1 (Napierian)
        return base + lmct


# --------------------------------------------------------------------------- #
# Onto the 480-band grid                                                        #
# --------------------------------------------------------------------------- #
def to_480(mac_fn, window_nm=(350.0, 800.0), uv_mode: str = "hold"):
    """Evaluate a MAC function on BioSNICAR's band centres with the window rules."""
    lo, hi = window_nm
    out = np.zeros(WVL_480_NM.size)
    inside = (WVL_480_NM >= lo) & (WVL_480_NM <= hi)
    out[inside] = mac_fn(WVL_480_NM[inside])
    below = WVL_480_NM < lo
    if uv_mode == "hold":
        out[below] = mac_fn(np.array([lo]))[0]
    elif uv_mode == "molecular":
        out[below] = mac_fn(WVL_480_NM[below])
    elif uv_mode != "zero":
        raise ValueError("uv_mode must be 'hold', 'molecular' or 'zero'")
    return np.clip(out, 0.0, None)


def water_k_480(biosnicar_root: str) -> np.ndarray:
    """Imaginary refractive index of water on the 480 grid (BioSNICAR's k file, 1-nm, 200-4999 nm)."""
    k = np.loadtxt(os.path.join(biosnicar_root, "data", "OP_data", "k_ice_480.csv"))
    return k[5::10]                 # same subsampling as biooptical_funcs.rescale_480band


# --------------------------------------------------------------------------- #
# Per-cell optical properties                                                   #
# --------------------------------------------------------------------------- #
@dataclass
class CellModel:
    geom: CellGeometry = field(default_factory=CellGeometry)
    c_internal_kg_m3: float = 50.0       # pigment mass / cell volume
    vacuole_fraction: float = 1.0
    water_volume_fraction: float = 0.625  # 0.59 * 1060/1000, as in BioSNICAR's bio-optical model
    n_real: float = 1.4                   # BioSNICAR default N_ALGAE (only used if g_mode='vd2014')
    g_mode: str = "fixed"                 # 'fixed' (g_fixed for every tier) or 'vd2014'
    g_fixed: float = 0.96                 # = BioSNICAR's empirical glacier-algae g (tier A)
    q_func: object = None                 # optional fast Q*(a, geom) (e.g. Phase 3 look-up table)
    vd_diagnostic: bool = True            # also compute the vd2014 SSA cross-check (slow, ~10 ms)

    @property
    def pigment_mass_per_cell_kg(self) -> float:
        return self.c_internal_kg_m3 * self.geom.volume_um3 * 1e-18

    def optics(self, mac480, k_water480, packaged: bool = True, scatter_from=None):
        """Return dict(ext_xsc, ss_alb, asm_prm, abs_xsc, q_star) on the 480 grid.

        packaged=False gives tier B (pigment absorbs as if dissolved); pass the
        packaged result as `scatter_from` so B keeps C's scattering cross-section.
        """
        from biosnicar.optical_properties.van_diedenhoven import calc_ssa_and_g

        g = self.geom
        V = g.volume_um3 * 1e-18                       # m^3
        A = g.projected_area_um2 * 1e-12               # m^2
        lam = WVL_480_UM * 1e-6                        # m
        m_pig = self.pigment_mass_per_cell_kg

        # intracellular absorption coefficients [m^-1]
        c_comp = self.c_internal_kg_m3 / self.vacuole_fraction
        a_pig = mac480 * c_comp
        a_wat = 4.0 * np.pi * k_water480 * self.water_volume_fraction / lam

        qf = self.q_func or q_star
        if self.vacuole_fraction >= 1.0:
            q = qf(a_pig + a_wat, g)
            abs_pig_pk = mac480 * m_pig * q
            abs_wat = a_wat * V * q
        else:
            q = qf(a_pig, g.scaled(self.vacuole_fraction))
            abs_pig_pk = mac480 * m_pig * q
            abs_wat = a_wat * V * qf(a_wat, g)

        abs_packaged = abs_pig_pk + abs_wat
        ext_geo = 2.0 * A
        if packaged:
            abs_x = abs_packaged
            sca = np.maximum(ext_geo - abs_packaged, 0.0)
        else:
            if scatter_from is None:
                raise ValueError("unpackaged tier needs the packaged optics for its scattering")
            sca = scatter_from["sca_xsc"]
            abs_x = mac480 * m_pig + abs_wat
        ext = sca + abs_x
        ssa = sca / ext

        # asymmetry parameter (van Diedenhoven 2014, as used by BioSNICAR for cells)
        a_tot = a_pig + a_wat
        k_cell = a_tot * lam / (4.0 * np.pi)
        ar = 1.0 if g.shape == "sphere" else (2.0 * g.radius) / g.length
        if self.vd_diagnostic or self.g_mode == "vd2014":
            ssa_vd, asym = calc_ssa_and_g(ar, g.volume_um3, g.projected_area_um2,
                                          np.full(WVL_480_UM.size, self.n_real), k_cell, WVL_480_UM)
        else:
            ssa_vd = asym = np.full(WVL_480_UM.size, np.nan)
        if self.g_mode == "fixed":
            # Cells embedded in ice/meltwater have a relative refractive index of only
            # ~1.05-1.10, so they are strongly forward scattering; the vd2014
            # parameterisation with n = 1.4 (relative to air) underestimates g.
            asym = np.full(WVL_480_UM.size, self.g_fixed)
        elif self.g_mode != "vd2014":
            raise ValueError("g_mode must be 'fixed' or 'vd2014'")
        return dict(ext_xsc=ext, ss_alb=np.clip(ssa, 1e-8, 1 - 1e-8), asm_prm=np.clip(asym, 0.0, 0.99),
                    abs_xsc=abs_x, sca_xsc=sca, q_star=q, ssa_vandiedenhoven=ssa_vd,
                    pigment_abs_xsc=(abs_pig_pk if packaged else mac480 * m_pig))


def model_a_optics(biosnicar_root: str, stem: str = "ice_algae_empirical_Chevrollier2023"):
    """Tier A: BioSNICAR's default empirical glacier-algae optical properties."""
    d = np.load(os.path.join(biosnicar_root, "data", "OP_data", "480band", "lap.npz"))
    ext, ssa, g = d[f"{stem}__ext_xsc"], d[f"{stem}__ss_alb"], d[f"{stem}__asm_prm"]
    return dict(ext_xsc=ext, ss_alb=ssa, asm_prm=g, abs_xsc=ext * (1.0 - ssa), sca_xsc=ext * ssa)


def describe(cell: CellModel) -> dict:
    d = asdict(cell)
    d.pop("geom")
    d.update(shape=cell.geom.shape, radius_um=cell.geom.radius, length_um=cell.geom.length,
             volume_um3=cell.geom.volume_um3, projected_area_um2=cell.geom.projected_area_um2,
             pigment_pg_per_cell=cell.pigment_mass_per_cell_kg * 1e15)
    return d
