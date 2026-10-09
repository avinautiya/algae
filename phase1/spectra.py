"""
Line-spectrum -> continuous spectrum -> Mass Absorption Coefficient (MAC).

Pure NumPy/SciPy/pandas: no quantum chemistry needed, so you can re-broaden
saved TD-DFT results (e.g. try FWHM = 0.25 / 0.35 eV) in seconds.

Physics
-------
For a transition i with oscillator strength f_i, the integrated molar
(decadic) absorption coefficient is (Hilborn, Am. J. Phys. 50, 982 (1982)):

    integral eps(nu~) d nu~  =  f_i / 4.319e-9      [L mol^-1 cm^-1 * cm^-1]

Changing variable to photon energy E in eV (nu~[cm^-1] = 8065.544 * E[eV]):

    integral eps(E) dE = f_i * K,   K = 1 / (4.319e-9 * 8065.544) = 2.8707e4 L mol^-1 cm^-1 eV

Each line is broadened with an area-normalised Gaussian in ENERGY space
(vibronic + inhomogeneous solvent broadening is approximately symmetric in
energy, not in wavelength):

    g_i(E) = 1/(sigma sqrt(2 pi)) * exp(-(E - E_i)^2 / (2 sigma^2)),
    sigma  = FWHM / (2 sqrt(2 ln 2))

    eps(E) = K * sum_i f_i g_i(E)                      [L mol^-1 cm^-1]

eps is a cross-section (not a density), so eps(lambda) = eps(E = hc/lambda)
with no Jacobian.  Conversion to mass absorption coefficient, with M in g/mol:

    1 L mol^-1 cm^-1 = 1e-3 m^3 mol^-1 * 1e2 m^-1 = 0.1 m^2 mol^-1
    MAC_decadic   [m^2 g^-1] = 0.1 * eps / M
    MAC_napierian [m^2 g^-1] = ln(10) * 0.1 * eps / M     (default)

The Napierian form is the one used by radiative-transfer snow/ice models
(e.g. BioSNICAR: absorption coefficient = MAC * mass concentration); set
napierian=False if you compare against decadic absorbance-based measurements.
Results are reported in m^2 kg^-1 (= 1000 * m^2 g^-1), BioSNICAR's unit.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

HC_EV_NM = 1239.841984              # h*c in eV*nm  (E[eV] = 1239.84 / lambda[nm])
EPS_INTEGRAL_K = 1.0 / (4.319e-9 * 8065.544)   # = 2.8707e4 L mol^-1 cm^-1 eV per unit f
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def nm_to_ev(lam_nm):
    return HC_EV_NM / np.asarray(lam_nm, dtype=float)


def ev_to_nm(e_ev):
    return HC_EV_NM / np.asarray(e_ev, dtype=float)


def gaussian_broaden(e_grid_ev, e_lines_ev, f_lines, fwhm_ev: float = 0.3):
    """Molar absorption coefficient eps(E) [L mol^-1 cm^-1] on `e_grid_ev`.

    Vectorised sum of area-normalised Gaussians (energy domain).
    """
    e_grid = np.asarray(e_grid_ev, dtype=float)[:, None]          # (n_grid, 1)
    e_i = np.asarray(e_lines_ev, dtype=float)[None, :]            # (1, n_states)
    f_i = np.asarray(f_lines, dtype=float)[None, :]
    sigma = fwhm_ev * FWHM_TO_SIGMA
    g = np.exp(-0.5 * ((e_grid - e_i) / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))
    return EPS_INTEGRAL_K * (f_i * g).sum(axis=1)


def epsilon_to_mac(eps_molar, molar_mass_g_mol: float, napierian: bool = True):
    """eps [L mol^-1 cm^-1] -> MAC [m^2 kg^-1]."""
    factor = np.log(10.0) if napierian else 1.0
    mac_m2_per_g = factor * 0.1 * np.asarray(eps_molar) / molar_mass_g_mol
    return 1000.0 * mac_m2_per_g


def fwhm_ev_to_nm(fwhm_ev: float, center_nm: float) -> float:
    """Approximate FWHM in nm of an energy-space Gaussian centred at center_nm."""
    e0 = HC_EV_NM / center_nm
    return float(ev_to_nm(e0 - fwhm_ev / 2) - ev_to_nm(e0 + fwhm_ev / 2))


# Wavelength window consumed downstream: phase2.tddft_calibration fits 260-750 nm extract MAC and
# normalises the HPLC shape and the Fe increment over 265-600 nm. Root coverage is judged here.
DOWNSTREAM_WINDOW_NM = (260.0, 750.0)
CAL_NORM_WINDOW_NM = (265.0, 600.0)


def root_count_sensitivity(energies_ev, osc, molar_mass_g_mol: float, fwhm_ev: float,
                           window_nm=DOWNSTREAM_WINDOW_NM, norm_nm=CAL_NORM_WINDOW_NM, drop: int = 5,
                           shift_ev: float = 0.0):
    """How much the broadened MAC in the downstream window changes when the highest `drop` roots are
    removed. A small change is evidence (not proof) that roots above the computed set would change the
    window little; a large change shows the window is NOT converged in the root count.

    Returns max relative pointwise change of MAC in window_nm (relative to the window maximum), the
    relative change of the integral over norm_nm, and the same at the window's short-wavelength edge."""
    e, f = np.asarray(energies_ev, float) + shift_ev, np.asarray(osc, float)
    lam = np.arange(window_nm[0], window_nm[1] + 0.5, 1.0)
    eg = nm_to_ev(lam)
    full = epsilon_to_mac(gaussian_broaden(eg, e, f, fwhm_ev), molar_mass_g_mol)
    n = e.size
    k = max(n - drop, 0)
    cut = epsilon_to_mac(gaussian_broaden(eg, e[:k], f[:k], fwhm_ev), molar_mass_g_mol) if k else 0 * full
    w = (lam >= norm_nm[0]) & (lam <= norm_nm[1])
    i_full, i_cut = np.trapezoid(full[w], lam[w]), np.trapezoid(cut[w], lam[w])
    return dict(n_roots=int(n), dropped=int(min(drop, n)), fwhm_ev=float(fwhm_ev), shift_ev=float(shift_ev),
                window_nm=list(window_nm), norm_window_nm=list(norm_nm),
                highest_root_ev=float(e.max()), edge_ev=float(HC_EV_NM / window_nm[0]),
                max_rel_change_in_window=float(np.max(np.abs(full - cut)) / max(full.max(), 1e-300)),
                rel_change_norm_integral=float(abs(i_full - i_cut) / max(i_full, 1e-300)),
                rel_change_at_short_edge=float(abs(full[0] - cut[0]) / max(full[0], 1e-300)))


def build_spectrum(energies_ev, osc_strengths, molar_mass_g_mol: float,
                   lam_min_nm: float = 300.0, lam_max_nm: float = 800.0,
                   step_nm: float = 1.0, fwhm_ev: float = 0.3, napierian: bool = True):
    """Discrete TD-DFT lines -> continuous spectrum + MAC.

    Returns
    -------
    lines_df : one row per excited state, columns
               [State, Wavelength_nm, Energy_eV, Oscillator_Strength, MAC_estimated]
               MAC_estimated = broadened MAC (m^2 kg^-1) evaluated at that state's
               wavelength (NaN if outside the grid range).
    spectrum_df : continuous spectrum on the wavelength grid, columns
               [Wavelength_nm, Energy_eV, Oscillator_Strength, Epsilon_L_mol-1_cm-1, MAC_estimated]
               Oscillator_Strength here is the broadened oscillator-strength density
               sum_i f_i g_i(E) in eV^-1 (dimensionless f per unit energy).
    arrays : dict of plain NumPy arrays (wavelength_nm, energy_ev, epsilon, mac_m2_kg,
             stick_wavelength_nm, stick_f) for plotting.
    """
    e_lines = np.asarray(energies_ev, dtype=float)
    f_lines = np.asarray(osc_strengths, dtype=float)

    lam = np.arange(lam_min_nm, lam_max_nm + 0.5 * step_nm, step_nm)
    e_grid = nm_to_ev(lam)
    eps = gaussian_broaden(e_grid, e_lines, f_lines, fwhm_ev)
    mac = epsilon_to_mac(eps, molar_mass_g_mol, napierian)
    f_density = eps / EPS_INTEGRAL_K

    spectrum_df = pd.DataFrame({
        "Wavelength_nm": lam,
        "Energy_eV": e_grid,
        "Oscillator_Strength": f_density,
        "Epsilon_L_mol-1_cm-1": eps,
        "MAC_estimated": mac,
    })

    lam_lines = ev_to_nm(e_lines)
    mac_at_lines = epsilon_to_mac(gaussian_broaden(e_lines, e_lines, f_lines, fwhm_ev),
                                  molar_mass_g_mol, napierian)
    in_range = (lam_lines >= lam_min_nm) & (lam_lines <= lam_max_nm)
    lines_df = pd.DataFrame({
        "State": np.arange(1, len(e_lines) + 1),
        "Wavelength_nm": lam_lines,
        "Energy_eV": e_lines,
        "Oscillator_Strength": f_lines,
        "MAC_estimated": np.where(in_range, mac_at_lines, np.nan),
    })

    arrays = dict(wavelength_nm=lam, energy_ev=e_grid, epsilon=eps, mac_m2_kg=mac,
                  stick_wavelength_nm=lam_lines, stick_f=f_lines)
    return lines_df, spectrum_df, arrays


def plot_spectrum(arrays, title: str = "", path: str | None = None, ax=None):
    """MAC(lambda) curve with oscillator-strength sticks on a twin axis."""
    import matplotlib
    if path is not None and ax is None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7, 4))
    lam = arrays["wavelength_nm"]
    ax.plot(lam, arrays["mac_m2_kg"], color="#7a2e0e", lw=2, label="MAC (broadened)")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel(r"MAC (m$^2$ kg$^{-1}$)")
    ax.set_xlim(lam.min(), lam.max())
    ax.set_ylim(bottom=0)

    ax2 = ax.twinx()
    sel = (arrays["stick_wavelength_nm"] >= lam.min()) & (arrays["stick_wavelength_nm"] <= lam.max())
    ax2.vlines(arrays["stick_wavelength_nm"][sel], 0, arrays["stick_f"][sel],
               color="#555555", lw=1.2, label="f (TD-DFT)")
    ax2.set_ylabel("Oscillator strength f")
    ax2.set_ylim(bottom=0, top=max(0.05, 1.15 * arrays["stick_f"][sel].max()) if sel.any() else 0.05)
    if title:
        ax.set_title(title)
    if own:
        fig.tight_layout()
        if path is not None:
            fig.savefig(path, dpi=600)
            plt.close(fig)
    return ax
