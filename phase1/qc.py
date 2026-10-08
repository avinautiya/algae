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


def optimize_geometry(mol, functional: str = "B3LYP", solvent: str | None = "pcm",
                      eps: float = WATER_EPS_STATIC, maxsteps: int = 100,
                      workdir: str = ".", tag: str = "opt", use_gpu: bool = False,
                      density_fit: bool = True, convergence: str = "gaussian", **mf_kwargs):
    """Ground-state geometry optimisation with geomeTRIC (TRIC internal coordinates).

    Every accepted step is appended to <workdir>/<tag>_traj.xyz and the latest
    geometry is written to <workdir>/<tag>_last.xyz, so a disconnected Colab
    session can be resumed by passing that file back in as the start geometry.

    convergence='gaussian' uses Gaussian's default thresholds
    (max force 4.5e-4, RMS force 3.0e-4, max step 1.8e-3, RMS step 1.2e-3 a.u.).
    Returns (optimised Mole, converged SCF object at that geometry, energy in Hartree,
             optimisation-converged flag).
    """
    from pyscf.geomopt.geometric_solver import kernel as geometric_kernel

    os.makedirs(workdir, exist_ok=True)
    traj = os.path.join(workdir, f"{tag}_traj.xyz")
    last = os.path.join(workdir, f"{tag}_last.xyz")
    open(traj, "w").close()
    t0 = time.time()

    mf = make_mf(mol, functional, solvent, eps, density_fit=density_fit, use_gpu=use_gpu, **mf_kwargs)

    def callback(envs):
        # geomeTRIC passes its local namespace; 'mol' is the current PySCF geometry.
        cur = envs["mol"]
        e = envs.get("energy", float("nan"))
        step = callback.n = getattr(callback, "n", 0) + 1
        xyz = _mole_to_xyz(cur, comment=f"step {step}  E = {e:.10f} Eh  "
                                         f"t = {time.time() - t0:.0f} s")
        with open(traj, "a") as fh:
            fh.write(xyz)
        with open(last, "w") as fh:
            fh.write(xyz)
        print(f"[{tag}] step {step:3d}  E = {e:.8f} Eh  elapsed {time.time() - t0:7.0f} s", flush=True)

    conv = {}
    if convergence == "gaussian":
        # gradients in Eh/Bohr, displacements in Angstrom (geomeTRIC units)
        conv = dict(convergence_energy=1e-6, convergence_grms=3.0e-4, convergence_gmax=4.5e-4,
                    convergence_drms=1.2e-3, convergence_dmax=1.8e-3)

    converged, mol_eq = geometric_kernel(mf, maxsteps=maxsteps, callback=callback, **conv)
    if not converged:
        print(f"WARNING [{tag}]: optimisation NOT converged in {maxsteps} steps; resume with "
              f"--start-xyz {last}", flush=True)

    # Final single point at the converged geometry (also gives MOs for TD-DFT).
    mf_eq = make_mf(mol_eq, functional, solvent, eps, density_fit=density_fit, use_gpu=use_gpu,
                    **mf_kwargs)
    e_eq = mf_eq.kernel()
    if not mf_eq.converged:
        raise RuntimeError("SCF at the optimised geometry did not converge")
    with open(os.path.join(workdir, f"{tag}_final.xyz"), "w") as fh:
        fh.write(_mole_to_xyz(mol_eq, comment=f"{functional}/{mol.basis} {solvent} E = {e_eq:.10f} Eh"))
    print(f"[{tag}] optimisation finished in {time.time() - t0:.0f} s, E = {e_eq:.10f} Eh", flush=True)
    return mol_eq, mf_eq, e_eq, converged


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
              conv_tol: float = 1e-6):
    """Vertical singlet->singlet excitations from a converged closed-shell SCF.

    For a closed-shell RKS reference PySCF's TDDFT/TDA solve only for singlet
    states by default (td.singlet = True).

    equilibrium_solvation=False (default): non-equilibrium linear-response PCM,
    the correct choice for vertical absorption.

    Returns dict with excitation energies (eV), wavelengths (nm), oscillator
    strengths (length gauge) and the TD object.
    """
    if not mf.converged:
        raise RuntimeError("SCF not converged; refusing to run TD-DFT on it")
    has_solvent = getattr(mf, "with_solvent", None) is not None
    if has_solvent:
        td = mf.TDA(equilibrium_solvation=equilibrium_solvation) if tda else \
            mf.TDDFT(equilibrium_solvation=equilibrium_solvation)
    else:
        td = mf.TDA() if tda else mf.TDDFT()
    if mol_spin(mf) == 0:
        td.singlet = True          # closed shell: singlet excitations only
    td.nstates = nstates
    td.conv_tol = conv_tol
    td.max_cycle = 200
    t0 = time.time()
    td.kernel()
    e = _to_numpy(td.e)
    f = _to_numpy(td.oscillator_strength(gauge="length"))
    conv = np.atleast_1d(_to_numpy(td.converged)).astype(bool)
    if not conv.all():
        print(f"WARNING: {np.sum(~conv)} of {len(conv)} TD-DFT roots not converged", flush=True)
    e_ev = e * HARTREE_TO_EV
    print(f"TD-DFT: {len(e_ev)} states, {e_ev.min():.3f}-{e_ev.max():.3f} eV, "
          f"{time.time() - t0:.0f} s", flush=True)
    return dict(energies_ev=e_ev, wavelengths_nm=1239.841984 / e_ev,
                osc_strengths=f, converged=conv, td=td)


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


def xyz_file_to_atom_block(path: str) -> str:
    """Read a (single- or multi-frame) .xyz file and return the LAST frame as an atom block."""
    with open(path) as fh:
        lines = [l.rstrip() for l in fh if l.strip()]
    frames, i = [], 0
    while i < len(lines):
        n = int(lines[i].split()[0])
        frames.append("\n".join(lines[i + 2:i + 2 + n]))
        i += 2 + n
    return frames[-1]


def set_threads(n: int | None = None):
    """Use all available cores (Colab: 2 vCPU on the free tier, 8 on High-RAM)."""
    n = n or os.cpu_count()
    lib.num_threads(n)
    return n
