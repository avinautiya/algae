"""
Content fingerprints, atomic writes and environment records for cached results.

A cached result is reused only if the fingerprint of everything it was computed from (inputs, settings,
the source code of the physics modules, the BioSNICAR revision and the numerical environment) matches
the fingerprint stored with it. A row count, a file name or a manually bumped version string is not
enough.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))

# modules whose code determines the optics / forward model; their content enters every physics fingerprint
PHYSICS_MODULES = ("phase1/spectra.py", "phase2/biosnicar_bridge.py", "phase2/cell_optics.py",
                   "phase2/pigment_packaging.py", "phase2/tddft_calibration.py", "phase2/empirical_data.py",
                   "phase3/forward_model.py", "phase4/emulator.py")


def _canon(obj):
    """JSON-serialisable canonical form; arrays by dtype, shape and a hash of their bytes."""
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {"__dc__": type(obj).__name__, **{f.name: _canon(getattr(obj, f.name)) for f in dataclasses.fields(obj)}}
    if isinstance(obj, np.ndarray):
        a = np.ascontiguousarray(obj)
        return {"__nd__": str(a.dtype), "shape": list(a.shape), "sha": hashlib.sha256(a.tobytes()).hexdigest()}
    if isinstance(obj, (np.floating, np.integer, np.bool_)):
        return obj.item()
    if isinstance(obj, dict):
        return {str(k): _canon(v) for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))}
    if isinstance(obj, (list, tuple)):
        return [_canon(v) for v in obj]
    if hasattr(obj, "to_numpy") and hasattr(obj, "columns"):           # DataFrame
        return {"__df__": list(map(str, obj.columns)), "values": _canon(obj.to_numpy())}
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return repr(obj)


def fingerprint(obj) -> str:
    return hashlib.sha256(json.dumps(_canon(obj), sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def code_fingerprint(paths=PHYSICS_MODULES) -> dict:
    out = {}
    for p in paths:
        full = os.path.join(ROOT, p)
        out[p] = file_sha256(full) if os.path.isfile(full) else None
    return out


def biosnicar_revision(root: str | None) -> str | None:
    if not root:
        return None
    try:
        return subprocess.check_output(["git", "-C", root, "rev-parse", "HEAD"], text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def environment() -> dict:
    import platform
    out = dict(python=platform.python_version(), numpy=np.__version__)
    for mod in ("scipy", "pandas", "emcee", "pyscf", "miepython", "rasterio"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:  # noqa: BLE001
            out[mod] = None
    return out


def phase1_fingerprint(phase1_dir: str | None) -> dict | None:
    """Content of a Phase 1 result directory that the optics consume (states, summary, marker)."""
    if not phase1_dir or not os.path.isdir(phase1_dir):
        return None
    files = sorted(f for f in os.listdir(phase1_dir)
                   if f.endswith(("_states.csv", "summary.json", "COMPLETE.json")))
    return {f: file_sha256(os.path.join(phase1_dir, f)) for f in files}


def empirical_data_fingerprint() -> str:
    d = os.path.join(ROOT, "data", "empirical")
    if not os.path.isdir(d):
        return "missing"
    h = hashlib.sha256()
    for f in sorted(os.listdir(d)):
        p = os.path.join(d, f)
        if os.path.isfile(p):
            h.update(f.encode())
            h.update(file_sha256(p).encode())
    return h.hexdigest()


def atomic_write_text(path: str, text: str):
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def atomic_to_csv(df, path: str, **kw):
    tmp = f"{path}.tmp.{os.getpid()}"
    df.to_csv(tmp, **kw)
    with open(tmp, "rb+") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def write_meta(path: str, fp: str, record: dict | None = None):
    """Sidecar <path>.meta.json: fingerprint of the inputs and the checksum of the result file."""
    meta = dict(fingerprint=fp, result_sha256=file_sha256(path), environment=environment(), record=record or {})
    atomic_write_text(path + ".meta.json", json.dumps(meta, indent=1, default=str))


def cached_ok(path: str, fp: str) -> tuple[bool, str]:
    """A cached result is valid only if its sidecar matches the expected fingerprint and the result file
    is the one the sidecar was written for."""
    side = path + ".meta.json"
    if not (os.path.isfile(path) and os.path.isfile(side)):
        return False, "no result or no provenance sidecar"
    try:
        meta = json.load(open(side))
    except (OSError, ValueError) as e:
        return False, f"sidecar unreadable: {e}"
    if meta.get("fingerprint") != fp:
        return False, "inputs changed (fingerprint mismatch)"
    if meta.get("result_sha256") != file_sha256(path):
        return False, "result file changed after it was written"
    return True, "ok"


if __name__ == "__main__":
    json.dump(dict(environment=environment(), code=code_fingerprint()), sys.stdout, indent=1)
