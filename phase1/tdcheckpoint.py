"""
Restartable linear-response (TDA / full TD) solves with validated checkpoints.

Design (see docs/repair_ledger.md, items P1-CKPT-*):

* The operator is fingerprinted: geometry (Bohr), atoms, charge/spin, basis and Cartesian convention,
  functional, RKS/UKS, density-fitting basis, integration grid, solvent model and response settings,
  TD vs TDA, spin channel, frozen orbitals, root count, software versions. The target residual
  tolerance is NOT part of the fingerprint (changing it permits reuse).
* The SCF orbitals used to define the MO-basis vectors are stored with the vectors. On resume the stored
  orbitals are compared with the current ones through the AO overlap; if the occupied and virtual
  subspaces agree (all singular values of the block overlaps > 1 - MO_SUBSPACE_TOL) the vectors are
  rotated into the current orbital basis (X' = U_oo^T X U_vv, likewise Y); otherwise the checkpoint
  is invalidated. Matching dimensions alone is never accepted.
* The installed PySCF (2.14) eigensolvers expose no iteration callback, so iteration-level progress is
  saved by running the solver in chunks of `chunk` Davidson iterations (td.max_cycle) and restarting
  from the Ritz vectors (a restarted Davidson; it converges to the same eigenpairs at the cost of extra
  iterations, measured in tests). This is chunk-level, not per-iteration, checkpointing.
* Three things are kept distinct: the latest restart vectors, the last stage at which ALL roots were
  genuinely converged (with its energies and oscillator strengths), and the final validated result.
* Writes are atomic (tmp + fsync + rename) with the previous valid checkpoint kept as `.prev`;
  a lock file prevents concurrent writers; a corrupt latest file falls back to `.prev`, and an
  unreadable pair is quarantined (renamed *.corrupt) and the solve restarts.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import time

import numpy as np

SCHEMA = 2
GEOM_TOL_BOHR = 1e-5          # a re-read .xyz can flip the 6th decimal of the rounded fingerprint geometry
MO_SUBSPACE_TOL = 1e-6
MO_ENERGY_TOL = 1e-6          # Eh: stored and current canonical orbital energies must agree
DEGENERATE_TOL = 1e-5         # Eh: orbitals closer than this may mix (degenerate rotations)


class CheckpointError(RuntimeError):
    pass


# --------------------------------------------------------------------------- fingerprint
def _round(a, d=6):
    return np.round(np.asarray(a, float), d).tolist()


def operator_fingerprint(mf, *, tda: bool, nstates: int, singlet: bool | None,
                         equilibrium_solvation: bool, frozen=None) -> dict:
    import pyscf
    mol = mf.mol
    ws = getattr(mf, "with_solvent", None)
    wdf = getattr(mf, "with_df", None)
    grids = getattr(mf, "grids", None)
    basis_repr = json.dumps(mol._basis, sort_keys=True, default=str)
    fp = dict(
        schema=SCHEMA,
        atoms=[mol.atom_symbol(i) for i in range(mol.natm)],
        coords_bohr=_round(mol.atom_coords(unit="Bohr"), 6),
        charge=int(mol.charge), spin=int(mol.spin),
        basis_sha=hashlib.sha256(basis_repr.encode()).hexdigest()[:16], cart=bool(mol.cart), nao=int(mol.nao),
        xc=str(getattr(mf, "xc", "HF")), scf_type=type(mf).__name__,
        restricted=not isinstance(mf.mo_occ, (tuple, list)) and np.asarray(mf.mo_occ).ndim == 1,
        df_auxbasis=None if wdf is None else str(getattr(wdf, "auxbasis", None)),
        grid_level=None if grids is None else int(getattr(grids, "level", -1)),
        grid_prune=None if grids is None else getattr(getattr(grids, "prune", None), "__name__", str(grids.prune)),
        solvent=None if ws is None else dict(type=type(ws).__name__, method=str(getattr(ws, "method", "")),
                                            eps=float(getattr(ws, "eps", np.nan)),
                                            lebedev_order=int(getattr(ws, "lebedev_order", -1))),
        equilibrium_solvation=bool(equilibrium_solvation),
        method="TDA" if tda else "RPA", singlet=singlet, frozen=frozen, nstates=int(nstates),
        pyscf=pyscf.__version__, numpy=".".join(np.__version__.split(".")[:2]))
    return fp


def fingerprint_hash(fp: dict) -> str:
    return hashlib.sha256(json.dumps(fp, sort_keys=True, default=str).encode()).hexdigest()


# --------------------------------------------------------------------------- orbitals and vectors
def equivalent_except_rounding(stored: dict, current: dict):
    """(True, max |dR| in Bohr) if two operator fingerprints differ ONLY through coordinate rounding:
    every non-geometry field identical, same atoms, max coordinate difference <= GEOM_TOL_BOHR. A geometry
    re-read from an .xyz file rounds differently from the in-memory optimised geometry (observed: one
    coordinate differing by 1e-6 Bohr), which would otherwise discard a valid checkpoint."""
    if not stored or set(stored) != set(current):
        return False, None
    if any(stored[k] != current[k] for k in stored if k != "coords_bohr"):
        return False, None
    a, b = np.asarray(stored["coords_bohr"], float), np.asarray(current["coords_bohr"], float)
    if a.shape != b.shape:
        return False, None
    dev = float(np.max(np.abs(a - b))) if a.size else 0.0
    return dev <= GEOM_TOL_BOHR, dev


def _spin_blocks(mf):
    """List of (C, occ, e) per spin channel."""
    occ = mf.mo_occ
    if isinstance(occ, (tuple, list)) or np.asarray(occ).ndim == 2:
        return [(np.asarray(mf.mo_coeff[s]), np.asarray(mf.mo_occ[s]), np.asarray(mf.mo_energy[s])) for s in (0, 1)]
    return [(np.asarray(mf.mo_coeff), np.asarray(mf.mo_occ), np.asarray(mf.mo_energy))]


def vectors_from_td(td, tda: bool):
    """Ritz vectors in the x0 format of td.kernel."""
    out = []
    for x, y in td.xy:
        xs = x if isinstance(x, tuple) else (x,)
        flat_x = np.hstack([np.ravel(np.asarray(v)) for v in xs])
        if tda:
            out.append(flat_x)
        else:
            ys = y if isinstance(y, tuple) else (y,)
            out.append(np.hstack([flat_x] + [np.ravel(np.asarray(v)) for v in ys]))
    return np.array(out)


def _block_shapes(mf):
    return [(int((occ > 0).sum()), int((occ == 0).sum())) for _, occ, _ in _spin_blocks(mf)]


def map_vectors(saved, mf, vecs, tda: bool):
    """Rotate MO-basis vectors from the saved orbitals into the current orbitals, or raise
    CheckpointError if the occupied/virtual subspaces differ."""
    S = mf.mol.intor_symmetric("int1e_ovlp")
    cur = _spin_blocks(mf)
    nspin = len(cur)
    if int(saved["nspin"]) != nspin:
        raise CheckpointError("spin channels differ")
    Us = []
    for s in range(nspin):
        C0, occ0, e0 = saved[f"mo_coeff_{s}"], saved[f"mo_occ_{s}"], saved[f"mo_energy_{s}"]
        C1, occ1, e1 = cur[s]
        if C0.shape != C1.shape or not np.array_equal(occ0 > 0, occ1 > 0):
            raise CheckpointError("orbital dimensions or occupations differ")
        if np.max(np.abs(np.asarray(e0) - np.asarray(e1))) > MO_ENERGY_TOL:
            raise CheckpointError("orbital energies differ (different SCF solution)")
        O = C0.T @ S @ C1
        # the response operator assumes canonical orbitals: only (near-)degenerate orbitals may mix
        nondeg = np.abs(np.asarray(e0)[:, None] - np.asarray(e1)[None, :]) > DEGENERATE_TOL
        if np.max(np.abs(O[nondeg]), initial=0.0) > 1e-4:
            raise CheckpointError("orbitals are not the same canonical set (non-degenerate mixing)")
        o, v = occ0 > 0, occ0 == 0
        Uoo, Uvv = O[np.ix_(o, o)], O[np.ix_(v, v)]
        for name, U in (("occupied", Uoo), ("virtual", Uvv)):
            sv = np.linalg.svd(U, compute_uv=False)
            if sv.min() < 1 - MO_SUBSPACE_TOL or sv.max() > 1 + MO_SUBSPACE_TOL:
                raise CheckpointError(f"{name} orbital subspace changed (singular values "
                                      f"{sv.min():.8f}-{sv.max():.8f})")
        Us.append((Uoo, Uvv))
    shapes = [(U[0].shape[0], U[1].shape[0]) for U in Us]
    size_x = sum(a * b for a, b in shapes)
    out = []
    for vec in vecs:
        parts = [vec[:size_x]] if tda else [vec[:size_x], vec[size_x:]]
        if not tda and parts[1].size != size_x:
            raise CheckpointError("vector length does not match the response space")
        if tda and vec.size != size_x:
            raise CheckpointError("vector length does not match the response space")
        new = []
        for part in parts:
            off = 0
            for (no, nv), (Uoo, Uvv) in zip(shapes, Us):
                blk = part[off:off + no * nv].reshape(no, nv)
                new.append((Uoo.T @ blk @ Uvv).ravel())
                off += no * nv
        out.append(np.hstack(new))
    return np.array(out)


# --------------------------------------------------------------------------- store
class CheckpointStore:
    """latest.npz (restart vectors + state), latest.prev.npz (previous valid), lock file."""

    def __init__(self, directory: str):
        self.dir = directory
        os.makedirs(directory, exist_ok=True)
        self.latest = os.path.join(directory, "latest.npz")
        self.prev = os.path.join(directory, "latest.prev.npz")
        self._lockf = None

    def __enter__(self):
        self._lockf = open(os.path.join(self.dir, ".lock"), "w")
        try:
            fcntl.flock(self._lockf, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise CheckpointError(f"checkpoint {self.dir} is locked by another process") from None
        self._lockf.write(str(os.getpid()))
        self._lockf.flush()
        return self

    def __exit__(self, *exc):
        fcntl.flock(self._lockf, fcntl.LOCK_UN)
        self._lockf.close()

    @staticmethod
    def _validate(d: dict, fp_hash: str):
        meta = json.loads(str(d["meta"]))
        if meta.get("schema") != SCHEMA:
            raise CheckpointError("schema mismatch")
        if meta.get("fp_hash") != fp_hash:
            raise CheckpointError("operator fingerprint mismatch")
        x = d["vectors"]
        if x.ndim != 2 or x.shape[0] != meta["nstates"] or not np.all(np.isfinite(x)):
            raise CheckpointError("restart vectors have the wrong shape or are not finite")
        for k in ("energies", "conv"):
            if d[k].shape != (meta["nstates"],):
                raise CheckpointError(f"{k} has the wrong shape")
        return meta

    def load(self, fp_hash: str, log=print):
        for path in (self.latest, self.prev):
            if not os.path.isfile(path):
                continue
            try:
                with np.load(path, allow_pickle=False) as z:
                    d = {k: z[k] for k in z.files}
                meta = self._validate(d, fp_hash)
                if path == self.prev:
                    log(f"[ckpt] latest checkpoint unusable; recovered previous checkpoint {path}")
                return d, meta
            except CheckpointError as e:
                if "fingerprint" in str(e):
                    log(f"[ckpt] {path}: {e}; checkpoint belongs to a different operator - not reused")
                    return None, None
                log(f"[ckpt] {path}: invalid ({e})")
            except Exception as e:  # noqa: BLE001 - corrupt / truncated file
                log(f"[ckpt] {path}: unreadable ({type(e).__name__}: {e})")
        for path in (self.latest, self.prev):          # quarantine unusable files
            if os.path.isfile(path):
                os.replace(path, path + f".corrupt.{int(time.time())}")
        return None, None

    def save(self, arrays: dict, meta: dict):
        tmp = os.path.join(self.dir, f".tmp.{os.getpid()}.npz")
        with open(tmp, "wb") as fh:
            np.savez(fh, meta=np.array(json.dumps(meta, default=str)), **arrays)
            fh.flush()
            os.fsync(fh.fileno())
        if os.path.isfile(self.latest):
            os.replace(self.latest, self.prev)
        os.replace(tmp, self.latest)
        dfd = os.open(self.dir, os.O_RDONLY)
        try:
            os.fsync(dfd)
        finally:
            os.close(dfd)


# --------------------------------------------------------------------------- solver
def default_stages(conv_tol):
    return [t for t in (1e-2, 1e-3, 1e-4) if t > conv_tol * (1 + 1e-9)] + [conv_tol]


def solve(td, mf, *, tda: bool, nstates: int, conv_tol: float, ckpt_dir: str, singlet=None,
          equilibrium_solvation=False, chunk: int = 5, max_total_cycles: int = 300,
          import_guess: str | None = None, log=print):
    """Run td to `conv_tol` in stages and chunks with checkpoints. Returns dict(status, energies_au,
    conv, stage_history, ...). status is 'converged' only if every root converged at conv_tol.
    import_guess: optional legacy .npz with an 'x0' array, used ONLY as an unvalidated initial guess
    (no stage credit) when no valid checkpoint exists."""
    fp = operator_fingerprint(mf, tda=tda, nstates=nstates, singlet=singlet,
                              equilibrium_solvation=equilibrium_solvation, frozen=getattr(td, "frozen", None))
    fph = fingerprint_hash(fp)
    stages = default_stages(conv_tol)
    with CheckpointStore(ckpt_dir) as store:
        # the stored fingerprint is read BEFORE it is overwritten, so a checkpoint whose operator differs only
        # by coordinate rounding can be recognised (restart vectors reused; stage credit is not)
        fpath = os.path.join(ckpt_dir, "fingerprint.json")
        stored = None
        if os.path.isfile(fpath):
            try:
                stored = json.load(open(fpath))
            except (OSError, ValueError):
                stored = None
        d, meta = store.load(fph, log)
        rounding_match = None
        if d is None and stored and stored.get("fp_hash") != fph:
            ok, dev = equivalent_except_rounding(stored.get("fp"), json.loads(json.dumps(fp, default=str)))
            if ok:
                d, meta = store.load(stored["fp_hash"], log)
                if d is not None:
                    rounding_match = dev
                    log(f"[ckpt] operator identical except coordinate rounding (max |dR| {dev:.1e} Bohr <= "
                        f"{GEOM_TOL_BOHR:g}); restart vectors reused, converged-stage credit NOT inherited")
        with open(fpath, "w") as fh:
            json.dump(dict(fp=fp, fp_hash=fph), fh, indent=1, default=str)
        x0, history, total, conv_done = None, [], 0, -1.0
        info = dict(resumed=False, invalidated=None, guess_imported=False, rounding_match_bohr=rounding_match,
                    resumed_as_guess=False)
        if d is not None:
            saved = dict(nspin=d["nspin"], **{k: d[k] for k in d if k.startswith("mo_")})
            try:
                x0 = list(map_vectors(saved, mf, d["vectors"], tda))
                history = meta.get("stage_history", [])
                total = int(meta.get("total_cycles", 0))
                conv_done = float(meta.get("converged_tol", -1.0)) if rounding_match is None else -1.0
                info["resumed"] = True
                log(f"[ckpt] resumed: {total} Davidson cycles done; all roots converged down to residual "
                    f"{conv_done if conv_done > 0 else 'none'}")
            except CheckpointError as e:
                if rounding_match is not None and np.asarray(d["vectors"]).shape[0] == nstates:
                    # operator identical except coordinate rounding, orbitals marginally different: the stored
                    # vectors are used ONLY as an unvalidated initial guess (as for import_guess) - the solve
                    # reconverges on the current operator to the full tolerance and earns no stage credit
                    x0 = list(np.asarray(d["vectors"], float))
                    history, total, conv_done = meta.get("stage_history", []), int(meta.get("total_cycles", 0)), -1.0
                    info.update(resumed=True, resumed_as_guess=True, invalidated=str(e))
                    log(f"[ckpt] orbitals differ marginally after coordinate rounding ({e}); stored vectors used as an "
                        f"UNVALIDATED initial guess only (no stage credit)")
                else:
                    log(f"[ckpt] stored orbitals incompatible with the current SCF ({e}); checkpoint invalidated")
                    x0, info["invalidated"] = None, str(e)
        if x0 is None and import_guess and os.path.isfile(import_guess):
            try:
                with np.load(import_guess) as z:
                    g = np.asarray(z["x0"], float)
                if g.shape[0] == nstates and np.all(np.isfinite(g)):
                    x0 = list(g)
                    info["guess_imported"] = True
                    log(f"[ckpt] using {import_guess} as an UNVALIDATED initial guess only (no stage credit)")
            except Exception as e:  # noqa: BLE001
                log(f"[ckpt] legacy guess {import_guess} unusable: {e}")
        mo = {}
        for s, (C, occ, e) in enumerate(_spin_blocks(mf)):
            mo.update({f"mo_coeff_{s}": C, f"mo_occ_{s}": occ, f"mo_energy_{s}": e})
        mo["nspin"] = np.array(len(_spin_blocks(mf)))

        def save(conv, energies):
            store.save(dict(vectors=np.asarray(x0), energies=np.asarray(energies, float),
                            conv=np.asarray(conv, bool), **mo),
                       dict(schema=SCHEMA, fp_hash=fph, nstates=nstates, total_cycles=total,
                            converged_tol=conv_done, stage_history=history, updated=time.time()))

        def start_vectors():
            """Restart vectors augmented with the solver's standard initial guess, so that a restart
            from Ritz vectors alone cannot lose a low-lying root (observed in tests without it)."""
            if x0 is None:
                return None
            guess = np.asarray(td.get_init_guess(mf, nstates))
            M = np.vstack([np.asarray(x0), guess])
            # orthonormal basis of the combined span; drop linearly dependent directions (small response
            # spaces, or guesses already contained in the Ritz space) which the eigensolvers reject
            _, sv, vt = np.linalg.svd(M, full_matrices=False)
            keep = sv > sv[0] * 1e-8
            if keep.sum() < nstates:
                return np.asarray(x0)
            return vt[keep]

        def run_chunk():
            """td.kernel with the augmented start set; if the eigensolver rejects it (singular trial
            space, seen for very small response spaces), retry with the plain Ritz vectors, then with a
            fresh standard guess. Each fallback is logged."""
            for label, xs in (("augmented", start_vectors), ("ritz", lambda: None if x0 is None else np.asarray(x0)),
                              ("fresh", lambda: None)):
                try:
                    td.kernel(x0=xs())
                    if label != "augmented":
                        log(f"[ckpt] eigensolver restarted with the {label} start set")
                    return
                except (RuntimeError, np.linalg.LinAlgError) as e:
                    log(f"[ckpt] eigensolver failed with the {label} start set ({e}); falling back")
            raise RuntimeError("eigensolver failed with every start set")

        status = "converged"
        if conv_done > 0 and conv_done <= conv_tol * (1 + 1e-9):
            # already complete: one verification pass on the current operator (also repopulates td)
            td.conv_tol, td.max_cycle = conv_tol, chunk
            run_chunk()
            conv = np.atleast_1d(np.asarray(td.converged, bool))
            log(f"[ckpt] verification pass: {conv.sum()}/{conv.size} roots converged at {conv_tol:g}")
            if not (conv.size == nstates and conv.all()):
                status = "verification_failed"
            return dict(status=status, fp=fp, fp_hash=fph, conv=conv, energies_au=np.asarray(td.e, float),
                        stage_history=history, total_cycles=total, stages=stages, **info)
        for tol in stages:
            if conv_done > 0 and conv_done <= tol * (1 + 1e-9):
                continue
            while True:
                td.conv_tol, td.max_cycle = tol, chunk
                t0 = time.time()
                run_chunk()
                total += chunk
                conv = np.atleast_1d(np.asarray(td.converged, bool))
                x0 = list(vectors_from_td(td, tda))
                energies = np.asarray(td.e, float)
                log(f"[ckpt] residual target {tol:g}: {conv.sum()}/{conv.size} roots converged after "
                    f"<= {total} cycles ({time.time() - t0:.0f} s this chunk)")
                if conv.size == nstates and conv.all():
                    conv_done = tol
                    f = np.asarray(td.oscillator_strength(gauge="length"), float)
                    history.append(dict(tol=tol, total_cycles=total, energies_au=energies.tolist(),
                                        osc=f.tolist()))
                    save(conv, energies)
                    break
                save(conv, energies)
                if total >= max_total_cycles:
                    status = "max_cycles_unconverged"
                    break
            if status != "converged":
                break
        return dict(status=status, fp=fp, fp_hash=fph, conv=conv, energies_au=np.asarray(td.e, float),
                    stage_history=history, total_cycles=total, stages=stages, **info)
