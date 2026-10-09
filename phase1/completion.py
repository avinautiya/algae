"""
Completion markers for Phase 1 runs (COMPLETE.json) and their validation.

A run directory is consumable downstream only through `validate_run`. Status levels:

  production            SCF converged; every requested root present and converged at the target
                        tolerance; finite positive excitation energies; finite non-negative oscillator
                        strengths; geometry status accepted (converged, or a measured and bounded
                        geometry sensitivity recorded); artifacts complete and checksummed.
  provisional_geometry  as production, but the geometry is not converged and its spectral sensitivity
                        has not been bounded yet. Usable only where a caller explicitly allows it, and
                        always reported.
  diagnostic_only       anything else (e.g. unconverged roots). Never consumed by calibration; no
                        COMPLETE.json is written, only DIAGNOSTIC.json.
"""

from __future__ import annotations

import hashlib
import json
import os

import numpy as np

MARKER, DIAG = "COMPLETE.json", "DIAGNOSTIC.json"
ACCEPTED_GEOMETRY = ("converged", "sensitivity_bounded")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def td_problems(energies_ev, osc, conv, nstates: int) -> list[str]:
    e, f, c = (np.asarray(x, float) for x in (energies_ev, osc, conv))
    probs = []
    if e.size != nstates:
        probs.append(f"{e.size} roots returned, {nstates} requested")
    if not np.all(np.asarray(conv, bool)):
        probs.append(f"{int(np.sum(~np.asarray(conv, bool)))} roots not converged")
    if not np.all(np.isfinite(e)) or np.any(e <= 0):
        probs.append("non-finite or non-positive excitation energies")
    if not np.all(np.isfinite(f)) or np.any(f < -1e-10):
        probs.append("non-finite or negative oscillator strengths")
    if e.size > 1 and np.any(np.diff(e) < -1e-8):
        probs.append("excitation energies not ascending")
    return probs


def geometry_status(start_xyz: str | None, optimised: bool | None, opt_converged: bool | None) -> dict:
    """Geometry acceptance. A start geometry is accepted only through a sidecar <xyz>.geometry.json
    (written by optimise runs or by an explicit sensitivity study)."""
    if optimised:
        return dict(status="converged" if opt_converged else "not_converged", source="this run")
    if start_xyz:
        side = start_xyz + ".geometry.json"
        if os.path.isfile(side):
            d = json.load(open(side))
            return dict(status=d.get("status", "unverified"), source=side, detail=d)
        return dict(status="unverified", source=start_xyz)
    return dict(status="unverified", source="built-in conformer (not optimised)")


def write_marker(outdir: str, summary: dict, artifacts: list[str], td_ok: bool, geom: dict):
    """Write COMPLETE.json (production / provisional_geometry) or DIAGNOSTIC.json, atomically."""
    if td_ok and geom["status"] in ACCEPTED_GEOMETRY:
        status, name = "production", MARKER
    elif td_ok:
        status, name = "provisional_geometry", MARKER
    else:
        status, name = "diagnostic_only", DIAG
    rec = dict(status=status, geometry=geom, settings={k: summary.get(k) for k in (
        "level", "basis", "cartesian_d", "solvent", "eps", "nstates", "tda", "fwhm_ev", "grid_level",
        "pcm_lebedev_order", "opt_functional", "td_conv_tol")},
        functionals=sorted(summary.get("tddft", {})),
        td=summary.get("tddft", {}), artifacts={os.path.basename(p): sha256(p) for p in artifacts})
    tmp = os.path.join(outdir, f".{name}.tmp")
    with open(tmp, "w") as fh:
        json.dump(rec, fh, indent=1, default=float)
    os.replace(tmp, os.path.join(outdir, name))
    other = os.path.join(outdir, DIAG if name == MARKER else MARKER)
    if os.path.exists(other):
        os.remove(other)
    return status


def validate_run(outdir: str, allow_provisional: bool = False, require: dict | None = None):
    """(ok, status, reasons). ok requires a marker whose artifact checksums match, production status
    (or provisional_geometry when allow_provisional), and matching settings for every key in require."""
    path = os.path.join(outdir, MARKER)
    if not os.path.isfile(path):
        st = "diagnostic_only" if os.path.isfile(os.path.join(outdir, DIAG)) else "missing"
        return False, st, [f"no {MARKER} in {outdir}"]
    try:
        rec = json.load(open(path))
    except (OSError, ValueError) as e:
        return False, "corrupt", [f"{MARKER} unreadable: {e}"]
    reasons = []
    for name, digest in rec.get("artifacts", {}).items():
        p = os.path.join(outdir, name)
        if not os.path.isfile(p):
            reasons.append(f"artifact {name} missing")
        elif sha256(p) != digest:
            reasons.append(f"artifact {name} changed after completion")
    for k, v in (require or {}).items():
        have = rec.get("settings", {}).get(k)
        if have != v:
            reasons.append(f"setting {k}={have!r}, required {v!r}")
    st = rec.get("status")
    if st == "provisional_geometry" and not allow_provisional:
        reasons.append("geometry not converged and its spectral sensitivity not yet bounded")
    elif st not in ("production", "provisional_geometry"):
        reasons.append(f"status {st}")
    return not reasons, st, reasons
