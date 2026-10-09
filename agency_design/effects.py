"""Conditional effect ranking from aligned joint ensemble draws.

This component neither retrieves algae nor generates scientifically valid draws.
The caller must supply independently checked, paired counterfactual model outputs.
"""
from __future__ import annotations

import numpy as np


def effects_from_albedo(without_algae, with_algae, flux_bin_w_m2):
    """Return signed W m-2 perturbations; layout (draw, region, wavelength).

    All three arrays must be aligned on draw/state/illumination and wavelength bins.
    Flux is BIN-INTEGRATED, not spectral density. No positivity clipping is used.
    Background dust/ice must be identical; the function cannot verify this metadata.
    """
    a0, a1, f = (np.asarray(v, float) for v in (without_algae, with_algae, flux_bin_w_m2))
    if a0.ndim != 3 or a0.shape != a1.shape or a0.shape != f.shape:
        raise ValueError("identical (draw, region, wavelength) arrays required")
    if not all(np.isfinite(v).all() for v in (a0, a1, f)):
        raise ValueError("nonfinite inputs; invalid pixels/draws must be resolved upstream")
    if np.any((a0 < 0) | (a0 > 1) | (a1 < 0) | (a1 > 1)) or np.any(f < 0):
        raise ValueError("albedo must be in [0,1] and downwelling flux nonnegative")
    return ((a0 - a1) * f).sum(axis=-1)


def aggregate_pixels(effect_draws, pixel_area_m2, *, quantity):
    """Preserve common draw IDs; no independent resampling of pixels.

    W m-2 becomes total W; J m-2 becomes total J. Totals and per-area effects must
    not be mixed in a ranking. Missing pixels are rejected rather than filled by zero.
    """
    q = np.asarray(effect_draws, float)
    area = np.asarray(pixel_area_m2, float)
    if quantity not in {"absorbed_sw_w_m2", "absorbed_energy_j_m2"}:
        raise ValueError("explicit supported effect quantity required")
    if q.ndim != 2 or area.shape != (q.shape[1],) or not np.isfinite(q).all():
        raise ValueError("effect shape (draw,pixel) and aligned area required")
    if not np.isfinite(area).all() or np.any(area <= 0):
        raise ValueError("finite positive pixel areas required")
    return dict(total=q @ area, mean_per_area=(q @ area) / area.sum(), area_m2=float(area.sum()),
                total_units="W" if quantity == "absorbed_sw_w_m2" else "J")


def rank_effects(effect_draws, region_ids, draw_ids, *, quantity, units,
                 minimum_effect, practical_difference, quality, probability_threshold=.95):
    """Rank conditional model effects and flag insufficient evidence.

    Layout (joint draw, region). minimum_effect/practical_difference are physical
    thresholds chosen BEFORE seeing maps, not learned cutoffs or species classes.
    Reported probabilities are conditional on supplied ensemble/model assumptions;
    they are not externally calibrated event probabilities or causal proof.
    """
    q = np.asarray(effect_draws, float)
    regions, draws = list(region_ids), list(draw_ids)
    supported = {"absorbed_sw": {"W m-2", "W"}, "absorbed_energy": {"J m-2", "J"}}
    if quantity not in supported or units not in supported[quantity]:
        raise ValueError("supported quantity/units required; this is not a melt-ranking module")
    if q.ndim != 2 or q.shape != (len(draws), len(regions)) or len(draws) < 2 or not regions:
        raise ValueError("at least two aligned draws and one region required")
    if len(set(regions)) != len(regions) or len(set(draws)) != len(draws):
        raise ValueError("unique region and joint-draw identifiers required")
    if not np.isfinite(q).all():
        raise ValueError("nonfinite draws must not disappear silently")
    if not np.isfinite([minimum_effect, practical_difference, probability_threshold]).all():
        raise ValueError("finite thresholds required")
    if minimum_effect < 0 or practical_difference < 0 or not .5 < probability_threshold < 1:
        raise ValueError("nonnegative physical thresholds and probability in (.5,1) required")
    required = {"observation_qc_pass", "in_domain", "retrieval_identifiable", "forward_validation_pass"}
    if set(quality) != set(regions):
        raise ValueError("quality metadata required for every region")
    rows = []
    for i, region in enumerate(regions):
        gates = quality[region]
        if not required <= gates.keys() or any(type(gates[k]) is not bool for k in required):
            raise ValueError("explicit boolean evidence gates required; missing gates cannot default to pass")
        failed = sorted(k for k in required if not gates[k])
        prob = float(np.mean(q[:, i] > minimum_effect))
        status = "unsupported" if failed else (
            "supported_model_effect" if prob >= probability_threshold else "effect_uncertain")
        rows.append(dict(region_id=region, median=float(np.median(q[:, i])),
                         interval95=np.quantile(q[:, i], [.025, .975]).tolist(),
                         p_exceeds_minimum_effect=prob, status=status, failed_gates=failed))
    pairs = []
    for i in range(len(regions)):
        for j in range(i + 1, len(regions)):
            delta = q[:, i] - q[:, j]
            pa, pb = float(np.mean(delta > practical_difference)), float(np.mean(-delta > practical_difference))
            eligible = rows[i]["status"] == rows[j]["status"] == "supported_model_effect"
            direction = regions[i] if pa >= probability_threshold else regions[j] if pb >= probability_threshold else None
            pairs.append(dict(a=regions[i], b=regions[j], difference_interval95=np.quantile(delta, [.025, .975]).tolist(),
                              p_a_exceeds_b_by_margin=pa, p_b_exceeds_a_by_margin=pb,
                              ranking="unsupported" if not eligible else "unresolved" if direction is None else "resolved",
                              higher_effect_region=direction if eligible else None))
    return dict(quantity=quantity, units=units, minimum_effect=minimum_effect,
                practical_difference=practical_difference, probability_threshold=probability_threshold,
                n_draws=len(draws), regions=rows, pairwise=pairs,
                interpretation="conditional ensemble effects; not calibrated causal/species/melt conclusions",
                total_order=None, warning="pairwise findings do not establish simultaneous global ranking confidence")
