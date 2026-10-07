"""
Empirical calibration of the Phase 1 TD-DFT pigment spectrum, and the empirical tier-D (Fe-complexed)
mass absorption coefficient.

Why: TD-DFT excitation energies, oscillator strengths and band widths carry method errors, and the
per-cell phenolic content we use (Williamson et al. 2020) is in PHENOL EQUIVALENTS (US EPA 420.1,
4-AAP), not in mass of the purpurogallin glucoside. Both are fixed by fitting the TD-DFT spectrum to
measured spectra of the same pigment:

  (1) shape: diode-array spectra of the chromatographically isolated (uncomplexed) glacier-algal
      phenolics, HPLC peaks 2-4 (Williamson et al. 2020 deposit), 265-600 nm, unit-area normalised;
  (2) magnitude: the mass absorption coefficient of whole phenolic extracts of S6 surface ice
      (53 samples, per kg phenol equivalents, with regression SE; Williamson et al. 2020), 260-750 nm,
      fitted in log space (so the visible, where sunlight peaks, weighs as much as the UV maximum),
      modelled as  E(l) = f * [ M(l; dE, w) + phi * I_M * D(l) ]
      where M is the TD-DFT MAC with every band shifted by dE (eV) and Gaussian FWHM w (eV), f absorbs
      the oscillator-strength error AND the glucoside -> phenol-equivalent mass conversion, phi in
      [0, 1] is the fraction of pigment complexed with Fe, D is the measured Fe-induced absorbance
      change of purpurogallin per unit integrated PG absorbance (Prochazkova et al. 2025, Fig. 4;
      empirical_data.fe_increment) and I_M the integral of M over the same 265-600 nm window, so the
      complex adds the same relative absorption to the glucoside chromophore as Fe adds to
      purpurogallin. Whole extracts contain the complexed pigment (dark), the isolated HPLC peak does
      not (yellowish) - Prochazkova et al. (2025).

Each data set gets a model-discrepancy SD fitted alongside (otherwise a stick spectrum that cannot
match every shoulder would give an over-confident posterior). The joint posterior of
(dE, w, f, beta, log s_shape, log s_mac) is sampled with emcee under flat priors on bounded ranges.

Outputs used downstream (all per kg of phenol equivalents, Napierian m^2 kg^-1):
  tier B/C  M_C(l) = f * M(l; dE, w)              calibrated uncomplexed pigment
  tier D    M_D(l) = f * [M(l; dE, w) + phi * I_M * D(l)]   pigment in its measured in-situ state
  Phase 3   marginal posteriors of dE, f, w, phi as the molecular-scale uncertainty PDFs.
"""

from __future__ import annotations

import dataclasses
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import empirical_data as ED  # noqa: E402

BOUNDS = dict(dE=(-1.0, 1.0), w=(0.08, 1.2), f=(1e-3, 1e3), phi=(0.0, 1.0), ls_shape=(-12.0, 0.0),
              ls_mac=(0.0, 14.0))
NAMES = list(BOUNDS)
MAC_RANGE_NM = (260.0, 750.0)
PHENOL_MOLAR_MASS = 94.11          # g/mol, C6H5OH: the standard of US EPA Method 420.1 (4-AAP)


def perturbed_mac(spec, dE: float = 0.0, f_scale: float = 1.0, fwhm: float | None = None):
    """MAC function of a Phase 1 spectrum with all bands shifted by dE (eV), oscillator strengths
    scaled by f_scale and Gaussian FWHM set to `fwhm` (eV). Continuous curves (demo input) are
    shifted, broadened in quadrature and scaled."""
    w = spec.fwhm_ev if fwhm is None else float(fwhm)
    if spec.has_sticks:
        s = dataclasses.replace(spec, energies_ev=spec.energies_ev + dE, osc=spec.osc * f_scale, fwhm_ev=w)
        return s.mac_at
    extra = float(np.sqrt(max(w ** 2 - spec.fwhm_ev ** 2, 0.0)))
    return lambda wl: f_scale * spec.mac_at(wl, shift_ev=dE, extra_fwhm_ev=extra)


@dataclasses.dataclass
class Calibration:
    spec: object
    samples: np.ndarray                 # (n, 6) posterior samples in NAMES order
    diagnostics: dict

    def mean(self, name):
        return float(self.samples[:, NAMES.index(name)].mean())

    def sd(self, name):
        return float(self.samples[:, NAMES.index(name)].std(ddof=1))

    def point(self):
        return {n: self.mean(n) for n in NAMES}

    def mac_C(self, wl, dE=None, f=None, w=None):
        p = self.point()
        dE, f, w = (p["dE"] if dE is None else dE), (p["f"] if f is None else f), (p["w"] if w is None else w)
        return perturbed_mac(self.spec, dE, f, w)(wl)

    def mac_D(self, wl, dE=None, f=None, w=None, phi=None):
        p = self.point()
        dE, f, w = (p["dE"] if dE is None else dE), (p["f"] if f is None else f), (p["w"] if w is None else w)
        phi = p["phi"] if phi is None else phi
        return complexed_mac(self.spec, dE, f, w, phi)(wl)

    def summary(self):
        out = {n: dict(mean=self.mean(n), sd=self.sd(n), q025=float(np.quantile(self.samples[:, i], 0.025)),
                       q975=float(np.quantile(self.samples[:, i], 0.975))) for i, n in enumerate(NAMES)}
        lf = np.log(self.samples[:, NAMES.index("f")])
        out["log_f"] = dict(mean=float(lf.mean()), sd=float(lf.std(ddof=1)))
        C = np.corrcoef(self.samples[:, :4].T)
        out["corr_dE_w_f_phi"] = np.round(C, 3).tolist()
        out.update(self.diagnostics)
        return out


_NORM_WL = np.arange(ED.FE_NORM_RANGE_NM[0], ED.FE_NORM_RANGE_NM[1] + 0.5, 1.0)


def complexed_mac(spec, dE, f, w, phi):
    """MAC of the pigment with a fraction phi complexed: f [M + phi I_M D] (see module docstring).
    Below 255 nm (no Fe data) the increment is held at its 255-nm value."""
    m = perturbed_mac(spec, dE, 1.0, w)
    i_m = np.trapezoid(m(_NORM_WL), _NORM_WL)
    d255 = float(ED.fe_increment([255.0])[0])

    def fn(wl):
        wl = np.asarray(wl, float)
        return f * (m(wl) + phi * i_m * np.nan_to_num(ED.fe_increment(wl), nan=d255))
    return fn


def _data():
    wl_h, S, sS = ED.chromophore_hplc_shape()
    wl_m, E, sE = ED.phenolic_extract_mac()
    keep = (wl_m >= MAC_RANGE_NM[0]) & (wl_m <= MAC_RANGE_NM[1])
    wl_m, E, sE = wl_m[keep], E[keep], sE[keep]
    return wl_h, S, sS, wl_m, E, sE, ED.fe_increment(wl_m)


def _shape_loglik(spec, wl_h, S, sS):
    def ll(th):
        dE, w, lsh = th
        Mh = perturbed_mac(spec, dE, 1.0, w)(wl_h)
        area = np.trapezoid(Mh, wl_h)
        if not np.isfinite(area) or area <= 0:
            return -np.inf
        vs = sS ** 2 + np.exp(2 * lsh)
        return -0.5 * np.sum((S - Mh / area) ** 2 / vs + np.log(2 * np.pi * vs))
    return ll


def _fit_magnitude(spec, dE, w, wl_m, E, sE, D):
    """Stage 2 for one (dE, w): maximum likelihood of (ln f, phi, ln s) with log-space residuals
    ln E - ln model, variance (sE/E)^2 + s^2; returns (theta, covariance of (ln f, phi))."""
    from scipy import optimize
    m1 = perturbed_mac(spec, dE, 1.0, w)
    M = m1(wl_m)
    i_m = np.trapezoid(m1(_NORM_WL), _NORM_WL)
    lnE, vr = np.log(E), (sE / E) ** 2

    def nll(t):
        mod = M + t[1] * i_m * D
        if np.any(mod <= 0):
            return 1e30
        v = vr + np.exp(2 * t[2])
        return 0.5 * np.sum((lnE - t[0] - np.log(mod)) ** 2 / v + np.log(v))

    best = None
    for phi0 in (0.0, 0.25, 0.5, 0.75, 1.0):
        lf0 = float(np.median(lnE - np.log(np.maximum(M + phi0 * i_m * D, 1e-30))))
        r = optimize.minimize(nll, [lf0, phi0, np.log(0.2)], method="L-BFGS-B",
                              bounds=[(np.log(BOUNDS["f"][0]), np.log(BOUNDS["f"][1])), BOUNDS["phi"], (-8, 3)])
        if best is None or r.fun < best.fun:
            best = r
    t0, h = best.x, np.array([1e-4, 1e-4, 1e-4])
    H = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            ea, eb = np.eye(3)[a] * h[a], np.eye(3)[b] * h[b]
            H[a, b] = (nll(t0 + ea + eb) - nll(t0 + ea - eb) - nll(t0 - ea + eb) + nll(t0 - ea - eb)) / (4 * h[a] * h[b])
    try:
        cov = np.linalg.inv(H)[:2, :2]
        if not np.all(np.isfinite(cov)) or np.any(np.diag(cov) < 0):
            raise np.linalg.LinAlgError
    except np.linalg.LinAlgError:
        cov = np.diag([1e-4, 1e-4])
    return t0, cov


def calibrate(spec, n_walkers=24, n_steps=2000, burn=800, n_stage2=200, draws_per=20, seed=0,
              verbose=True) -> Calibration:
    """Modular ('cut') posterior. Stage 1: (dE, w) from the isolated-chromophore shape alone (emcee).
    Stage 2: for draws of stage 1, (f, phi) from the extract MAC (log-space maximum likelihood with a
    fitted discrepancy SD, Gaussian approximation; phi kept in [0, 1]). Band position and width are
    properties of the uncomplexed chromophore, so the extract (which also holds Fe complexes and other
    phenolics) is not allowed to move them."""
    import emcee
    import logging
    logging.getLogger("emcee").setLevel(logging.ERROR)
    from scipy import optimize
    wl_h, S, sS, wl_m, E, sE, D = _data()
    ll = _shape_loglik(spec, wl_h, S, sS)
    lo = np.array([BOUNDS["dE"][0], BOUNDS["w"][0], BOUNDS["ls_shape"][0]])
    hi = np.array([BOUNDS["dE"][1], BOUNDS["w"][1], BOUNDS["ls_shape"][1]])
    lp = lambda t: ll(t) if np.all(t > lo) and np.all(t < hi) else -np.inf  # noqa: E731
    best = max(((lp(np.array([d, w, np.log(0.1 * S.max())])), d, w) for d in np.arange(-0.9, 0.91, 0.05)
                for w in (0.12, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0)), key=lambda x: x[0])
    x0 = optimize.minimize(lambda t: -lp(t) if np.isfinite(lp(t)) else 1e300,
                           [best[1], best[2], np.log(0.1 * S.max())], method="Nelder-Mead").x
    rng = np.random.default_rng(seed)
    p0 = np.clip(x0 + np.array([0.01, 0.01, 0.05]) * rng.normal(size=(n_walkers, 3)), lo + 1e-6, hi - 1e-6)
    smp = emcee.EnsembleSampler(n_walkers, 3, lp)
    smp.run_mcmc(p0, n_steps, progress=False)
    ch1 = smp.get_chain(discard=burn, flat=True)
    try:
        tau = smp.get_autocorr_time(discard=burn, quiet=True).tolist()
    except Exception:  # noqa: BLE001
        tau = None
    pick = ch1[rng.choice(len(ch1), size=min(n_stage2, len(ch1)), replace=False)]
    rows, lsm = [], []
    for dE, w, lsh in pick:
        t, cov = _fit_magnitude(spec, dE, w, wl_m, E, sE, D)
        z = rng.multivariate_normal(t[:2], cov, size=draws_per)
        z[:, 1] = np.clip(z[:, 1], 0.0, 1.0)
        for lf, ph in z:
            rows.append([dE, w, np.exp(lf), ph, lsh, t[2]])
        lsm.append(t[2])
    ch = np.array(rows)
    pm = ch.mean(axis=0)
    Mh = perturbed_mac(spec, pm[0], 1.0, pm[1])(wl_h)
    Sh = Mh / np.trapezoid(Mh, wl_h)
    Em = complexed_mac(spec, pm[0], pm[2], pm[1], pm[3])(wl_m)
    Ec = complexed_mac(spec, pm[0], pm[2], pm[1], 0.0)(wl_m)
    r2 = lambda o, m: float(1 - np.sum((o - m) ** 2) / np.sum((o - o.mean()) ** 2))  # noqa: E731
    vis = (wl_m >= 400) & (wl_m <= 700)
    lo6 = np.array([BOUNDS[n][0] for n in NAMES])
    hi6 = np.array([BOUNDS[n][1] for n in NAMES])
    diag = dict(spectrum=getattr(spec, "source", ""), method="cut posterior (shape -> dE, w; extract -> f, phi)",
                shape_r2=r2(S, Sh), mac_r2_log=r2(np.log(E), np.log(Em)),
                mac_r2_log_400_700=r2(np.log(E[vis]), np.log(Em[vis])),
                mac_rms_rel_400_700=float(np.sqrt(np.mean((Em[vis] / E[vis] - 1) ** 2))),
                acceptance_stage1=float(np.mean(smp.acceptance_fraction)), autocorr_steps_stage1=tau,
                at_bound=[NAMES[i] for i in range(3)
                          if min(pm[i] - lo6[i], hi6[i] - pm[i]) < 0.02 * (hi6[i] - lo6[i])],
                fe_share_of_absorption_400_700=float(
                    1 - np.trapezoid(Ec[vis], wl_m[vis]) / np.trapezoid(Em[vis], wl_m[vis])))
    # stoichiometric reading of f: if the 4-AAP assay counts one glucoside molecule as one phenol, the
    # mass conversion kg glucoside per kg phenol equivalent is M_spec / M_phenol, and the remaining
    # factor f / (M_spec / M_phenol) is the TD-DFT oscillator-strength (intensity) error alone.
    diag["phenol_molar_mass"] = PHENOL_MOLAR_MASS
    diag["f_stoichiometric_1to1"] = float(spec.molar_mass / PHENOL_MOLAR_MASS)
    diag["f_over_stoichiometric"] = float(pm[2] / (spec.molar_mass / PHENOL_MOLAR_MASS))
    cal = Calibration(spec, ch, diag)
    if verbose:
        sm = cal.summary()
        print(f"TD-DFT calibration ({diag['spectrum']}): dE = {sm['dE']['mean']:+.3f} +/- {sm['dE']['sd']:.3f} eV, "
              f"FWHM = {sm['w']['mean']:.3f} +/- {sm['w']['sd']:.3f} eV, f = {sm['f']['mean']:.3g} "
              f"+/- {sm['f']['sd']:.2g}, complexed fraction phi = {sm['phi']['mean']:.3f} +/- {sm['phi']['sd']:.3f}; "
              f"shape R2 {diag['shape_r2']:.3f}, extract ln-MAC R2 {diag['mac_r2_log']:.3f} "
              f"(400-700 nm rms rel. error {diag['mac_rms_rel_400_700']:.2f}); f / (M/M_phenol) = "
              f"{diag['f_over_stoichiometric']:.2f}")
        if diag["at_bound"]:
            print(f"  WARNING: posterior mean at a bound for {diag['at_bound']} - the TD-DFT spectrum does not "
                  "resemble the measured pigment spectrum; check the Phase 1 run")
    return cal


def calibrate_core(spec_l1, verbose=True):
    """Diagnostic: shift/width of the Level 1 (purpurogallin core) spectrum fitted to the measured
    absorbance of purpurogallin in water (Prochazkova et al. 2025, Fig. 4), 265-600 nm, by least
    squares on unit-area shapes. Returns dict(dE, w, r2)."""
    from scipy import optimize
    w_all, a_pg, _ = ED.fe_purpurogallin_absorbance()
    k = (w_all >= 265) & (w_all <= 600) & np.isfinite(a_pg)
    wl, A = w_all[k], np.clip(a_pg[k], 0, None)
    A = A / np.trapezoid(A, wl)

    def res(t):
        M = perturbed_mac(spec_l1, t[0], 1.0, t[1])(wl)
        a = np.trapezoid(M, wl)
        return (A - M / a) if a > 0 else np.full_like(A, 1e3)

    best = min(((np.sum(res([d, w]) ** 2), d, w) for d in np.arange(-0.9, 0.91, 0.05) for w in (0.2, 0.35, 0.5, 0.7)))
    r = optimize.least_squares(res, [best[1], best[2]], bounds=([-1.0, 0.08], [1.0, 1.2]))
    r2 = float(1 - np.sum(r.fun ** 2) / np.sum((A - A.mean()) ** 2))
    if verbose:
        print(f"Level 1 core vs measured purpurogallin (water): dE = {r.x[0]:+.3f} eV, FWHM = {r.x[1]:.3f} eV, "
              f"shape R2 = {r2:.3f}")
    return dict(dE=float(r.x[0]), w=float(r.x[1]), r2=r2)


def plot_calibration(cal: Calibration, figdir: str, F):
    """Fig. S0: (a) HPLC-isolated chromophore shape vs calibrated TD-DFT shape;
    (b) S6 extract MAC vs the calibrated model with and without the Fe-complexed fraction."""
    import matplotlib.pyplot as plt
    wl_h, S, sS, wl_m, E, sE, _ = _data()
    p = cal.point()
    Mh = perturbed_mac(cal.spec, p["dE"], 1.0, p["w"])(wl_h)
    Mr = perturbed_mac(cal.spec, 0.0, 1.0, None)(wl_h)
    fig, ax = plt.subplots(1, 2, figsize=(F.DOUBLE_COL, 2.7), constrained_layout=True)
    ax[0].fill_between(wl_h, S - 2 * sS, S + 2 * sS, color=F.INK_2, alpha=0.2, lw=0)
    ax[0].plot(wl_h, S, color=F.INK, label="HPLC peaks 2-4 (measured)")
    ax[0].plot(wl_h, Mr / np.trapezoid(Mr, wl_h), ls=":", color="#d55e00", label="TD-DFT, raw")
    ax[0].plot(wl_h, Mh / np.trapezoid(Mh, wl_h), color="#d55e00", label="TD-DFT, calibrated")
    ax[0].set_xlabel("Wavelength (nm)")
    ax[0].set_ylabel(r"Normalised absorbance (nm$^{-1}$)")
    ax[0].set_title(rf"a  shape: $\Delta E$={p['dE']:+.2f} eV, FWHM={p['w']:.2f} eV", loc="left")
    ax[0].legend(loc="upper right")
    ax[1].fill_between(wl_m, (E - 2 * sE) / 1e3, (E + 2 * sE) / 1e3, color=F.INK_2, alpha=0.2, lw=0)
    ax[1].plot(wl_m, E / 1e3, color=F.INK, label="S6 extract (measured)")
    ax[1].plot(wl_m, cal.mac_C(wl_m) / 1e3, color="#0072b2", label="tier C: calibrated, uncomplexed")
    ax[1].plot(wl_m, cal.mac_D(wl_m) / 1e3, ls="--", color="#009e73",
               label=rf"tier D: + Fe complex, $\phi$={p['phi']:.2f}")
    ax[1].set_yscale("log")
    ax[1].set_ylim(max(np.nanmin(E[E > 0]) / 1e3 / 3, 1e-3), None)
    ax[1].set_xlabel("Wavelength (nm)")
    ax[1].set_ylabel(r"MAC (10$^3$ m$^2$ kg$^{-1}$ phenol eq.)")
    ax[1].set_title(rf"b  magnitude: f={p['f']:.3g}", loc="left")
    ax[1].legend(loc="lower left")
    F.save(fig, figdir, "FigS0_tddft_calibration")


_CACHE = {}


def cached_calibration(spec, **kw) -> Calibration:
    """calibrate() once per spectrum per process (keyed by source and stick data)."""
    key = (getattr(spec, "source", ""), spec.name,
           None if spec.energies_ev is None else tuple(np.round(spec.energies_ev, 6)),
           None if spec.osc is None else tuple(np.round(spec.osc, 8)))
    if key not in _CACHE:
        _CACHE[key] = calibrate(spec, **kw)
    return _CACHE[key]
