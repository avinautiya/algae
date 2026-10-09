"""Versioned, non-pickle NPZ spectral interface with JSON provenance sidecar."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from .coupling import spectral_budget, _array

SCHEMA_VERSION = "algae-glacier-spectra-1"


def export_spectra(path, wavelength_nm, albedo_without_algae, albedo_with_algae,
                   flux_bin_w_m2, *, provenance):
    """Time/state/member axes remain explicit in provenance; last axis is wavelength.

    Export evaluated member spectra, not a claimed universal LUT. The host must
    match illumination, ice/dust and state coordinates; extrapolation is unsupported.
    Existing outputs are not overwritten. Supply draw IDs to preserve joint uncertainty.
    """
    dst = Path(path)
    if dst.suffix != ".npz":
        raise ValueError("export path must end in .npz")
    side = dst.with_suffix(".json")
    if dst.exists() or side.exists():
        raise FileExistsError("choose a new output path; exports are immutable")
    required = {"source_commit", "input_sha256", "biosnicar_commit", "state_coordinates",
                "leading_dimensions", "illumination", "calibration_draw_ids", "validation_status"}
    if not required <= provenance.keys():
        raise ValueError(f"provenance must include {sorted(required)}")
    w = _array(wavelength_nm, "wavelength_nm", 0)
    a0 = np.asarray(albedo_without_algae, dtype=float)
    a1 = np.asarray(albedo_with_algae, dtype=float)
    f = np.asarray(flux_bin_w_m2, dtype=float)
    spectral_budget(a0, a1, f)
    if w.ndim != 1 or len(w) != a0.shape[-1] or np.any(np.diff(w) <= 0):
        raise ValueError("wavelengths must be strictly increasing and match the spectral axis")
    if len(provenance["leading_dimensions"]) != a0.ndim - 1:
        raise ValueError("name every leading time/state/member dimension")
    metadata = dict(schema_version=SCHEMA_VERSION, provenance=provenance,
                    units={"wavelength_nm": "nm", "albedo": "1", "flux_bin_w_m2": "W m-2 per bin"},
                    shape=list(a0.shape), interpolation="none; evaluated states only")
    # Validate JSON before writing any artifact.
    json.dumps(metadata, allow_nan=False)
    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("xb") as handle:
        np.savez_compressed(handle, wavelength_nm=w, albedo_without_algae=a0,
                            albedo_with_algae=a1, flux_bin_w_m2=f)
    metadata["artifact_sha256"] = hashlib.sha256(dst.read_bytes()).hexdigest()
    with side.open("x") as handle:
        handle.write(json.dumps(metadata, indent=2, allow_nan=False) + "\n")
    return metadata


def load_spectra(path):
    dst = Path(path)
    metadata = json.loads(dst.with_suffix(".json").read_text())
    if metadata["schema_version"] != SCHEMA_VERSION:
        raise ValueError("unsupported schema")
    if hashlib.sha256(dst.read_bytes()).hexdigest() != metadata["artifact_sha256"]:
        raise ValueError("spectral export checksum mismatch")
    with np.load(dst, allow_pickle=False) as archive:
        arrays = {k: archive[k].copy() for k in archive.files}
    spectral_budget(arrays["albedo_without_algae"], arrays["albedo_with_algae"], arrays["flux_bin_w_m2"])
    return arrays, metadata
