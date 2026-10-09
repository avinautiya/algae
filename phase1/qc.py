"""
PySCF quantum-chemistry layer for Phase 1:

  * restricted Kohn-Sham ground state (B3LYP or CAM-B3LYP / 6-31G(d))
  * implicit water: IEF-PCM (eps = 78.4) or SMD (solvent = 'water')
  * geometry optimisation with geomeTRIC (checkpointed every step)
  * TD-DFT (full linear response, RPA) with non-equilibrium PCM for vertical
    absorption: the fast (electronic) solvent response uses the optical
    dielectric constant eps_inf = n^2 = 1.78 of water.

Everything optionally runs on a GPU through gpu4pyscf (`use_gpu=True`), which is
what makes the 48-atom Level 2 glycoside feasible inside a Colab session.
"""

from __future__ import annotations

import os
import time

import numpy as np
from pyscf import dft, lib

HARTREE_TO_EV = 27.211386245988

# Functional keywords.  'B3LYPG' is B3LYP with VWN(RPA) correlation, i.e. exactly
# the Gaussian 'B3LYP' definition (PySCF's own 'B3LYP' uses VWN5 and differs by
# a few mHartree).  Use 'B3LYP' instead if you benchmark against ORCA/Turbomole.
FUNCTIONALS = {
    "B3LYP": "B3LYPG",
    "CAM-B3LYP": "CAMB3LYP",
}

WATER_EPS_STATIC = 78.4      # task specification (78.3553 in Gaussian/PySCF tables)


def make_mf(mol, functional: str = "B3LYP", solvent: str | None = "pcm",
            eps: float = WATER_EPS_STATIC, density_fit: bool = True,
            grids_level: int = 3, conv_tol: float = 1e-9, use_gpu: bool = False,
            chkfile: str | None = None, lebedev_order: int = 29):
    """Build a restricted Kohn-Sham object with implicit solvation.

    functional : 'B3LYP' or 'CAM-B3LYP' (or any raw libxc/PySCF xc string)
    solvent    : 'pcm' -> IEF-PCM with dielectric `eps`
                 'smd' -> SMD with solvent='water' (eps taken from the SMD table, 78.3553)
                 None  -> gas phase
    density_fit: RI-J/K with the def2-universal-jkfit auxiliary basis; ~3-5x faster,
                 error in excitation energies < 0.01 eV.
    """
    xc = FUNCTIONALS.get(functional.upper(), functional)
    # closed shell -> restricted KS; open shell (e.g. high-spin Fe(III), Level 3) -> unrestricted KS
    mf = dft.RKS(mol, xc=xc) if mol.spin == 0 else dft.UKS(mol, xc=xc)
    if density_fit:
        mf = mf.density_fit(auxbasis="def2-universal-jkfit")
    if use_gpu:
        # Same construction order as gpu4pyscf's own drivers: DF -> to_gpu -> solvent.
        # Requires `pip install gpu4pyscf-cuda12x` on a CUDA runtime.
        mf = mf.to_gpu()
    mf.grids.level = grids_level          # level 3 ~ (75,302) pruned, Gaussian "FineGrid"-like
    mf.conv_tol = conv_tol
    mf.max_cycle = 200
    mf.init_guess = "minao"
    mf.chkfile = chkfile
    if chkfile and os.path.isfile(chkfile):
        # restart (container restarts): start from the last saved orbitals, projected onto this geometry
        mf.init_guess = "chkfile"

    if solvent is not None:
        solvent = solvent.lower()
        if solvent == "pcm":
            mf = mf.PCM()
            mf.with_solvent.method = "IEF-PCM"
            mf.with_solvent.eps = eps
        elif solvent == "smd":
            mf = mf.SMD()
            mf.with_solvent.solvent = "water"
        else:
            raise ValueError(f"unknown solvent model {solvent!r}")
        mf.with_solvent.lebedev_order = lebedev_order   # 29 -> 302 points/sphere (PySCF default)
    return mf


GAUSSIAN_CRITERIA = dict(convergence_energy=1e-6, convergence_grms=3.0e-4, convergence_gmax=4.5e-4,
                         convergence_drms=1.2e-3, convergence_dmax=1.8e-3)
CRITERIA_UNITS = dict(convergence_energy="Eh (energy change)", convergence_grms="Eh/Bohr",
                      convergence_gmax="Eh/Bohr", convergence_drms="Angstrom", convergence_dmax="Angstrom")


def optimize_geometry(mol, functional: str = "B3LYP", solvent: str | None = "pcm",
                      eps: float = WATER_EPS_STATIC, maxsteps: int = 100,
                      workdir: str = ".", tag: str = "opt", use_gpu: bool = False,
                      density_fit: bool = True, convergence: str = "gaussian", **mf_kwargs):
    """Ground-state geometry optimisation with geomeTRIC (TRIC internal coordinates).

    History is preserved across restarts: every gradient evaluation is APPENDED to
    <tag>_traj.xyz and to <tag>_steps.jsonl (energy, energy change, RMS/max Cartesian gradient in
    Eh/Bohr, RMS/max displacement from the previous evaluation in Angstrom, segment id). The latest
    geometry is <tag>_last.xyz. The optimiser's own convergence decision and every criterion are
    written to <tag>_opt_record.json.

    <tag>_final.xyz (the validated optimised geometry) is written ONLY if geomeTRIC reports
    convergence; otherwise the endpoint is written to <tag>_unconverged_endpoint.xyz and the record
    says so. (geomeTRIC's internal Hessian is not restorable through PySCF, so a restart begins a new
    segment with a fresh Hessian; the segments are recorded.)

    convergence='gaussian' uses Gaussian's default thresholds (GAUSSIAN_CRITERIA).
    Returns (Mole at the final/endpoint geometry, converged SCF there, energy in Eh, converged flag).
    """
    import json
    from pyscf.geomopt.geometric_solver import kernel as geometric_kernel

    os.makedirs(workdir, exist_ok=True)
    traj = os.path.join(workdir, f"{tag}_traj.xyz")
    last = os.path.join(workdir, f"{tag}_last.xyz")
    steps = os.path.join(workdir, f"{tag}_steps.jsonl")
    record_path = os.path.join(workdir, f"{tag}_opt_record.json")
    segment = 1
    if os.path.isfile(steps):
        with open(steps) as fh:
            segment = 1 + max((json.loads(l).get("segment", 0) for l in fh if l.strip()), default=0)
    t0 = time.time()
    mf = make_mf(mol, functional, solvent, eps, density_fit=density_fit, use_gpu=use_gpu, **mf_kwargs)
    state = dict(n=0, prev_xyz=None, prev_e=None, rows=[])

    def callback(envs):
        cur = envs["mol"]
        e = float(envs.get("energy", float("nan")))
        g = np.asarray(envs.get("gradients", np.full((cur.natm, 3), np.nan)), float)
        xyz_a = cur.atom_coords(unit="Angstrom")
        state["n"] += 1
        row = dict(segment=segment, step=state["n"], energy_Eh=e,
                   dE_Eh=None if state["prev_e"] is None else e - state["prev_e"],
                   grms=float(np.sqrt(np.mean(g ** 2))), gmax=float(np.max(np.abs(g))),
                   drms=None if state["prev_xyz"] is None else float(np.sqrt(np.mean((xyz_a - state["prev_xyz"]) ** 2))),
                   dmax=None if state["prev_xyz"] is None else float(np.max(np.abs(xyz_a - state["prev_xyz"]))),
                   elapsed_s=round(time.time() - t0, 1))
        state["prev_xyz"], state["prev_e"] = xyz_a, e
        state["rows"].append(row)
        frame = _mole_to_xyz(cur, comment=f"segment {segment} step {state['n']}  E = {e:.10f} Eh")
        with open(traj, "a") as fh:
            fh.write(frame)
        with open(steps, "a") as fh:
            fh.write(json.dumps(row) + "\n")
        _atomic_write(last, frame)
        print(f"[{tag}] seg {segment} step {state['n']:3d}  E = {e:.8f}  grms {row['grms']:.2e} "
              f"gmax {row['gmax']:.2e}  elapsed {time.time() - t0:7.0f} s", flush=True)

    conv = dict(GAUSSIAN_CRITERIA) if convergence == "gaussian" else {}
    converged, mol_eq = geometric_kernel(mf, maxsteps=maxsteps, callback=callback, **conv)
    mf_eq = make_mf(mol_eq, functional, solvent, eps, density_fit=density_fit, use_gpu=use_gpu, **mf_kwargs)
    e_eq = mf_eq.kernel()
    if not mf_eq.converged:
        raise RuntimeError("SCF at the optimised geometry did not converge")
    last_row = state["rows"][-1] if state["rows"] else {}
    rec = dict(tag=tag, functional=functional, basis=str(mol.basis), solvent=solvent, segment=segment,
               optimizer="geomeTRIC (via pyscf.geomopt)", converged=bool(converged), maxsteps=maxsteps,
               thresholds=conv, units=CRITERIA_UNITS, evaluations_this_segment=state["n"],
               last_evaluation=last_row, final_energy_Eh=float(e_eq),
               note="geomeTRIC's convergence flag is authoritative; the per-step values are recorded "
                    "for audit (displacements here are between successive gradient evaluations).")
    _atomic_write(record_path, json.dumps(rec, indent=1, default=float))
    xyz_out = os.path.join(workdir, f"{tag}_final.xyz" if converged else f"{tag}_unconverged_endpoint.xyz")
    _atomic_write(xyz_out, _mole_to_xyz(mol_eq, comment=f"{functional}/{mol.basis} {solvent} E = {e_eq:.10f} Eh "
                                                         f"converged={bool(converged)}"))
    if not converged:
        print(f"WARNING [{tag}]: optimisation NOT converged in {maxsteps} steps; endpoint -> {xyz_out}; "
              f"resume with --start-xyz {last}", flush=True)
    print(f"[{tag}] optimisation segment finished in {time.time() - t0:.0f} s, E = {e_eq:.10f} Eh", flush=True)
    return mol_eq, mf_eq, e_eq, converged


def _atomic_write(path: str, text: str):
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def xtb_preoptimize(mol, solvent: str | None = "water", gtol: float = 2e-4, maxiter: int = 2000):
    """Cheap GFN2-xTB / ALPB(water) pre-optimisation (seconds), via `pip install tblite`.

    Moves the MMFF starting structure close to the DFT minimum (planarises the
    benzotropolone, settles the O-H...O=C hydrogen-bond network), which typically
    halves the number of expensive DFT optimisation steps. Returns a new Mole.
    """
    from scipy.optimize import minimize
    from tblite.interface import Calculator

    numbers = mol.atom_charges()
    x0 = mol.atom_coords(unit="Bohr").ravel()

    def relax(method):
        calc = Calculator(method, numbers, x0.reshape(-1, 3), charge=float(mol.charge), uhf=mol.spin)
        calc.set("verbosity", 0)
        calc.set("max-iter", 500)
        if solvent:
            calc.add("alpb-solvation", solvent)

        def fun(x):
            calc.update(x.reshape(-1, 3))
            res = calc.singlepoint()
            return res.get("energy"), res.get("gradient").ravel()

        o = minimize(fun, x0, jac=True, method="L-BFGS-B", options=dict(gtol=gtol, maxiter=maxiter, maxcor=50))
        return o, np.abs(fun(o.x)[1]).max()

    t0 = time.time()
    method = "GFN2-xTB"
    try:
        opt, gmax = relax(method)
    except Exception as e:  # noqa: BLE001 - e.g. high-spin Fe(III) dications: GFN2 SCF may not converge
        print(f"[xtb] GFN2-xTB failed ({e}); falling back to GFN1-xTB", flush=True)
        method = "GFN1-xTB"
        opt, gmax = relax(method)
    xtb_preoptimize.method = method
    print(f"[xtb] {method}/ALPB({solvent}) pre-opt: E = {opt.fun:.6f} Eh, "
          f"max|g| = {gmax:.1e}, {opt.nit} iterations, {time.time() - t0:.0f} s", flush=True)
    return mol.set_geom_(opt.x.reshape(-1, 3), unit="Bohr", inplace=False)


def harmonic_check(mf, workdir: str = ".", tag: str = "freq"):
    """Optional: analytic PCM Hessian -> harmonic frequencies (cm^-1).

    Confirms the optimised structure is a true minimum (no imaginary modes).
    Expensive: ~3-6x an optimisation step per atom-batch; skip in tight sessions.
    """
    from pyscf.hessian import thermo

    h = mf.Hessian().kernel()
    res = thermo.harmonic_analysis(mf.mol, h)
    freqs = np.asarray(res["freq_wavenumber"])
    np.savetxt(os.path.join(workdir, f"{tag}_frequencies_cm-1.txt"), np.real_if_close(freqs))
    n_imag = int(np.sum(np.iscomplex(freqs) & (np.abs(np.imag(freqs)) > 1e-6)))
    print(f"[{tag}] {n_imag} imaginary frequencies", flush=True)
    return freqs, n_imag


def run_tddft(mf, nstates: int = 30, tda: bool = False, equilibrium_solvation: bool = False,
              conv_tol: float = 1e-5, ckpt_dir: str | None = None, chunk: int = 5,
              max_total_cycles: int = 400, import_guess: str | None = None):
    """Vertical excitations from a converged SCF (closed shell: singlets only).

    conv_tol is the Davidson residual-norm tolerance. It does not by itself guarantee a stated
    excitation-energy accuracy; the stage history (energies and oscillator strengths at each converged
    residual level) is returned so the effect of the tolerance can be measured.

    With ckpt_dir the solve is restartable and validated (phase1/tdcheckpoint.py); status is
    'converged' only if every root converged at conv_tol. Without ckpt_dir a single plain solve is run.

    equilibrium_solvation=False (default): non-equilibrium linear-response PCM (vertical absorption).

    Returns dict(energies_ev, wavelengths_nm, osc_strengths, converged, status, td, ...). Callers must
    check `status`/`converged` (see completion.td_problems) before using the numbers in production.
    """
    if not mf.converged:
        raise RuntimeError("SCF not converged; refusing to run TD-DFT on it")
    has_solvent = getattr(mf, "with_solvent", None) is not None
    if has_solvent:
        td = mf.TDA(equilibrium_solvation=equilibrium_solvation) if tda else \
            mf.TDDFT(equilibrium_solvation=equilibrium_solvation)
    else:
        td = mf.TDA() if tda else mf.TDDFT()
    singlet = None
    if mol_spin(mf) == 0:
        td.singlet = singlet = True          # closed shell: singlet excitations only
    td.nstates = nstates
    td.conv_tol = conv_tol
    td.max_cycle = 200
    t0 = time.time()
    extra = {}
    if ckpt_dir:
        import tdcheckpoint
        r = tdcheckpoint.solve(td, mf, tda=tda, nstates=nstates, conv_tol=conv_tol, ckpt_dir=ckpt_dir,
                               singlet=singlet, equilibrium_solvation=equilibrium_solvation, chunk=chunk,
                               max_total_cycles=max_total_cycles, import_guess=import_guess)
        status = r["status"]
        extra = dict(fingerprint=r["fp"], fp_hash=r["fp_hash"], stage_history=r["stage_history"],
                     total_cycles=r["total_cycles"], resumed=r["resumed"], invalidated=r["invalidated"])
    else:
        td.kernel()
        status = "converged" if np.all(np.atleast_1d(_to_numpy(td.converged))) else "unconverged"
    e = _to_numpy(td.e)
    f = _to_numpy(td.oscillator_strength(gauge="length"))
    conv = np.atleast_1d(_to_numpy(td.converged)).astype(bool)
    if not conv.all():
        print(f"WARNING: {np.sum(~conv)} of {len(conv)} TD-DFT roots not converged - results are "
              f"diagnostic only", flush=True)
    e_ev = e * HARTREE_TO_EV
    print(f"TD-DFT ({'TDA' if tda else 'RPA'}): {len(e_ev)} states, {e_ev.min():.3f}-{e_ev.max():.3f} eV, "
          f"status {status}, {time.time() - t0:.0f} s", flush=True)
    return dict(energies_ev=e_ev, wavelengths_nm=1239.841984 / e_ev,
                osc_strengths=f, converged=conv, status=status, td=td, **extra)


def mol_spin(mf):
    return getattr(getattr(mf, "mol", None), "spin", 0)


def dominant_transitions(td, state: int, threshold: float = 0.1):
    """Leading occupied->virtual contributions (|coefficient|^2 * 2 > threshold).
    For unrestricted (open-shell) references the alpha and beta channels are listed separately
    (labels a/b) with weights x^2 - y^2, which then sum to ~1 over both spins."""
    x, y = td.xy[state]
    if isinstance(x, (tuple, list)):                      # UKS: ((xa, xb), (ya, yb))
        out = []
        for spin, xs, ys in zip("ab", x, y if isinstance(y, (tuple, list)) else (0.0, 0.0)):
            xs = _to_numpy(xs)
            ys = 0.0 if np.isscalar(ys) else _to_numpy(ys)
            w = xs ** 2 - (ys ** 2 if not np.isscalar(ys) else 0.0)
            nocc = xs.shape[0]
            for i, a in zip(*np.where(np.abs(w) > threshold / 2)):
                hole = f"HOMO{spin}" + (f"-{nocc - 1 - i}" if i != nocc - 1 else "")
                part = f"LUMO{spin}" + (f"+{a}" if a else "")
                out.append((hole, part, float(w[i, a])))
        return sorted(out, key=lambda t: -t[2])
    x = _to_numpy(x)
    y = _to_numpy(y) if not np.isscalar(y) else 0.0
    w = 2.0 * (x ** 2 - (y ** 2 if not np.isscalar(y) else 0.0))   # weights sum to ~1
    nocc = x.shape[0]
    out = []
    for i, a in zip(*np.where(np.abs(w) > threshold)):
        hole = "HOMO" + (f"-{nocc - 1 - i}" if i != nocc - 1 else "")
        part = "LUMO" + (f"+{a}" if a else "")
        out.append((hole, part, float(w[i, a])))
    return sorted(out, key=lambda t: -t[2])


def _to_numpy(a):
    return a.get() if hasattr(a, "get") else np.asarray(a)


def _mole_to_xyz(mol, comment: str = "") -> str:
    coords = mol.atom_coords(unit="Angstrom")
    lines = [str(mol.natm), comment]
    for i in range(mol.natm):
        x, y, z = coords[i]
        lines.append(f"{mol.atom_pure_symbol(i):2s} {x:14.8f} {y:14.8f} {z:14.8f}")
    return "\n".join(lines) + "\n"


class XYZError(ValueError):
    pass


def read_xyz_frames(path: str, recover: bool = True):
    """Parse a single- or multi-frame XYZ file. Each frame: atom count, ONE comment line (which may be
    blank), then that many 'symbol x y z' lines. Returns (frames, warnings) with frames a list of
    (comment, [(symbol, x, y, z), ...]). A malformed or truncated frame after at least one valid frame
    is dropped with a warning if recover=True (last-valid-frame recovery); otherwise XYZError."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    frames, warnings, i = [], [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1                       # blank separator lines between frames / trailing newline
            continue
        try:
            n = int(lines[i].split()[0])
            if n <= 0:
                raise ValueError
        except (ValueError, IndexError):
            msg = f"line {i + 1}: expected an atom count, got {lines[i]!r}"
            if frames and recover:
                warnings.append(msg + "; trailing content ignored")
                break
            raise XYZError(msg) from None
        body = lines[i + 2:i + 2 + n]
        atoms = []
        try:
            if i + 1 >= len(lines) or len(body) < n:
                raise XYZError(f"frame at line {i + 1} truncated ({len(body)} of {n} atom lines)")
            for k, l in enumerate(body):
                parts = l.split()
                if len(parts) < 4:
                    raise XYZError(f"line {i + 3 + k}: malformed coordinate line {l!r}")
                x, y, z = (float(v) for v in parts[1:4])
                if not all(np.isfinite([x, y, z])):
                    raise XYZError(f"line {i + 3 + k}: non-finite coordinate")
                atoms.append((parts[0], x, y, z))
        except (XYZError, ValueError) as e:
            if frames and recover:
                warnings.append(f"{e}; recovered the last valid frame")
                break
            raise XYZError(str(e)) from None
        frames.append((lines[i + 1], atoms))
        i += 2 + n
    if not frames:
        raise XYZError(f"{path}: no frames")
    return frames, warnings


def xyz_file_to_atom_block(path: str, frame: int = -1) -> str:
    """Atom block of one frame (default: the last valid frame) of a .xyz file."""
    frames, warnings = read_xyz_frames(path)
    for w in warnings:
        print(f"WARNING [{path}]: {w}", flush=True)
    return "\n".join(f"{s} {x:.10f} {y:.10f} {z:.10f}" for s, x, y, z in frames[frame][1])


def set_threads(n: int | None = None):
    """Use all available cores (Colab: 2 vCPU on the free tier, 8 on High-RAM)."""
    n = n or os.cpu_count()
    lib.num_threads(n)
    return n
