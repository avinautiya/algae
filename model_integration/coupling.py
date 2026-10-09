"""Conservative spectral coupling and explicitly conditional surface-melt diagnostics.

Spectral flux inputs are BIN-INTEGRATED W m-2, not W m-2 nm-1.
Arrays have (..., wavelength) layout; time integration uses explicit interval durations.
Paired spectra must retain identical ice, dust, illumination and geometry.
"""
from __future__ import annotations

import numpy as np

LATENT_HEAT_J_KG = 334000.0


def _array(value, name, low=None, high=None):
    a = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be finite; missing data cannot be silently filled")
    if low is not None and np.any(a < low):
        raise ValueError(f"{name} must be >= {low}")
    if high is not None and np.any(a > high):
        raise ValueError(f"{name} must be <= {high}")
    return a


def spectral_budget(albedo_without_algae, albedo_with_algae, flux_bin_w_m2):
    """Absorbed SW and signed algae perturbation (last axis integrated).

    No positivity clipping: scattering can brighten some ice configurations.
    This is column absorption, NOT a guarantee that all heat melts the surface.
    """
    a0 = _array(albedo_without_algae, "albedo_without_algae", 0, 1)
    a1 = _array(albedo_with_algae, "albedo_with_algae", 0, 1)
    f = _array(flux_bin_w_m2, "flux_bin_w_m2", 0)
    if a0.shape != a1.shape or a0.shape != f.shape or a0.ndim < 1:
        raise ValueError("paired albedos and bin flux must have identical (..., wavelength) shapes")
    sw = f.sum(axis=-1)
    q0, q1 = ((1 - a) * f for a in (a0, a1))
    return dict(sw_down_w_m2=sw, absorbed_without_w_m2=q0.sum(axis=-1),
                absorbed_with_w_m2=q1.sum(axis=-1),
                delta_absorbed_w_m2=((a0 - a1) * f).sum(axis=-1))


def couple_host_albedo(host_albedo, modeled_without, modeled_with, *, mode,
                       host_algae_treatment):
    """Explicit replacement or additive anomaly; never infer the host's algae treatment.

    Replacement swaps the FULL albedo parameterization and requires matched nonalgal
    conditions. Anomaly correction is allowed only on an algae-free host. No clipping
    is used to hide an incompatible baseline. Validation of matching conditions remains
    the caller's responsibility.
    """
    h = _array(host_albedo, "host_albedo", 0, 1)
    a0 = _array(modeled_without, "modeled_without", 0, 1)
    a1 = _array(modeled_with, "modeled_with", 0, 1)
    if h.shape != a0.shape or h.shape != a1.shape:
        raise ValueError("host and paired modeled albedos must have identical shapes")
    if host_algae_treatment not in {"absent", "explicit", "implicit_observed"}:
        raise ValueError("declare host_algae_treatment: absent, explicit, or implicit_observed")
    if mode == "replace":
        return a1.copy()
    if mode != "anomaly":
        raise ValueError("mode must be replace or anomaly")
    if host_algae_treatment != "absent":
        raise ValueError("adding algae to an explicit/observed-albedo host would double count")
    return _array(h + a1 - a0, "corrected host albedo", 0, 1)


def mix_subpixels(albedo, area_fraction):
    """Area mix already-solved albedos, never average abundance before nonlinear RT.

    Shape (patch, ..., wavelength); shared irradiance assumed. If irradiance varies
    between patches, area-average their absorbed fluxes instead.
    """
    a = _array(albedo, "albedo", 0, 1)
    w = _array(area_fraction, "area_fraction", 0, 1)
    if a.ndim < 2 or w.ndim != 1 or len(w) != a.shape[0]:
        raise ValueError("one area fraction is required per patch")
    if not np.isclose(w.sum(), 1, rtol=0, atol=1e-10):
        raise ValueError("area fractions must sum to one")
    return np.tensordot(w, a, axes=(0, 0))


def integrate_intervals(flux_w_m2, dt_seconds):
    """Interval-MEAN flux -> J m-2. Axis 0 is time; gaps must be resolved upstream."""
    q = _array(flux_w_m2, "flux_w_m2")
    dt = _array(dt_seconds, "dt_seconds", 0)
    if q.ndim < 1 or dt.ndim != 1 or len(dt) != q.shape[0] or np.any(dt == 0):
        raise ValueError("positive explicit duration required for each time interval")
    return np.tensordot(dt, q, axes=(0, 0))


def melt_diagnostic(absorbed_sw_w_m2, non_sw_net_w_m2, dt_seconds, *,
                    initial_cold_content_j_m2, surface_at_melting_point):
    """Offline energy-limited potential melt, kg m-2 (= mm water equivalent).

    All inputs shape (time,). Positive non-SW is toward surface (LW, sensible,
    latent, rain and conductive boundary terms, excluding SW). Negative net energy
    replenishes the cold-content reservoir. No liquid-water retention, refreezing,
    runoff, penetrative heating or temperature solver is included. A host SEB must
    provide the melting-point eligibility flag. This is not independently validated
    actual melt and must not be labeled runoff or surface mass balance.
    """
    sw = _array(absorbed_sw_w_m2, "absorbed_sw_w_m2", 0)
    other = _array(non_sw_net_w_m2, "non_sw_net_w_m2")
    dt = _array(dt_seconds, "dt_seconds", 0)
    eligible = np.asarray(surface_at_melting_point)
    cold = _array(initial_cold_content_j_m2, "initial_cold_content_j_m2", 0)
    if sw.ndim != 1 or other.shape != sw.shape or dt.shape != sw.shape or np.any(dt <= 0):
        raise ValueError("SW, non-SW and positive durations must be identical 1D time arrays")
    if eligible.dtype != bool or eligible.shape != sw.shape or cold.ndim != 0:
        raise ValueError("eligibility must be a boolean time array; initial cold content a scalar")
    reservoir = float(cold)
    melt, withheld, history = [], [], []
    for energy, can_melt in zip((sw + other) * dt, eligible):
        if energy <= 0:
            reservoir -= energy
            remaining = 0.0
        else:
            warming = min(reservoir, energy)
            reservoir -= warming
            remaining = energy - warming
        melt.append(remaining / LATENT_HEAT_J_KG if can_melt else 0.0)
        withheld.append(0.0 if can_melt else remaining)
        history.append(reservoir)
    return dict(potential_melt_mm_we=np.array(melt), cold_content_j_m2=np.array(history),
                unallocated_energy_j_m2=np.array(withheld))


def paired_melt_diagnostic(without_sw, with_sw, non_sw, durations, **kwargs):
    """Same meteorology/initial state, separate evolving cold reservoirs for each run."""
    base = melt_diagnostic(without_sw, non_sw, durations, **kwargs)
    algae = melt_diagnostic(with_sw, non_sw, durations, **kwargs)
    return dict(without=base, with_algae=algae,
                delta_potential_melt_mm_we=algae["potential_melt_mm_we"] - base["potential_melt_mm_we"])
