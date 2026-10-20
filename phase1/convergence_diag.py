"""Bounded convergence diagnostics for the restartable TD-DFT solves (P9-CONV).

Read-only with respect to production results: checkpoints are COPIED into a diagnostic directory before any
solve; production SCF chk files are copied before an SCF restart (PySCF writes to its chkfile).

Subcommands
  config     numerical configuration actually in force (installed PySCF solver, search-space sizes, limits)
  residuals  per-root energies and TRUE residual norms of saved checkpoint vectors on the current operator,
             plus state tracking between checkpoints by RPA-metric amplitude overlaps and cluster subspace angles
  segment    run a bounded solve segment from a copy of a checkpoint through the production code path
             (tdcheckpoint.solve) with a given chunk length; records per-cycle per-root residuals (parsed from
             the solver's DEBUG1 log), operator applications (vectors passed through the response operator)
             and wall time per operator call

Residual convention: identical to pyscf.tdscf._lr_eig.real_eig (the solver in use): for a Ritz pair
(omega, X, Y) normalised to X.X - Y.Y = 1,
  r = || (A X + B Y - omega X ,  B X + A Y + omega Y) ||_2 .
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import molecules  # noqa: E402
import qc  # noqa: E402
import tdcheckpoint as TC  # noqa: E402

V2 = os.path.join(HERE, "results", "v2")
DIAG = os.path.join(HERE, "results", "convergence_diag")
SYSTEMS = {
    # the stalled production system (L1_FULL) and a small representative variant (same molecule, functional,
    # solvent and response type; smaller basis) for implementation tests
    "L1": dict(level="level1", xyz=os.path.join(V2, "L1_FULL", "opt_B3LYP_final.xyz"), basis="6-31g*",
               chk_dir=os.path.join(V2, "L1_FULL"), max_memory=2100),
    "L1_small": dict(level="level1", xyz=os.path.join(V2, "L1_FULL", "opt_B3LYP_final.xyz"), basis="3-21g",
                     chk_dir=None, max_memory=2100),
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def build(system: str, func: str, workdir: str, verbose=0, max_memory=None):
    """mol + converged mf exactly as run_phase1 builds them (DF def2-universal-jkfit, grid level 3, Lebedev 29,
    IEF-PCM eps 78.4, Cartesian 6-31G*), restarting the SCF from a COPY of the production chk file."""
    s = SYSTEMS[system]
    os.makedirs(workdir, exist_ok=True)
    mol = molecules.get_mole(s["level"], basis=s["basis"], xyz=qc.xyz_file_to_atom_block(s["xyz"]), cart=True,
                             max_memory=max_memory or s["max_memory"], verbose=verbose)
    chk = os.path.join(workdir, f"scf_{func}.chk")
    if s["chk_dir"] and not os.path.isfile(chk):
        src = os.path.join(s["chk_dir"], f"scf_{func}.chk")
        if os.path.isfile(src):
            shutil.copy2(src, chk)
    mf = qc.make_mf(mol, func, "pcm", qc.WATER_EPS_STATIC, density_fit=True, chkfile=chk, grids_level=3,
                    lebedev_order=29)
    t0 = time.time()
    mf.kernel()
    if not mf.converged:
        raise RuntimeError("SCF not converged")
    # orbital gradient with the SOLVATED Fock matrix (PySCF's default get_grad omits the PCM potential for
    # _Solvation objects: hcore + get_veff without v_solvent)
    g = mf.get_grad(mf.mo_coeff, mf.mo_occ, mf.get_fock(dm=mf.make_rdm1()))
    scf = dict(energy=float(mf.e_tot), seconds=time.time() - t0, conv_tol=mf.conv_tol,
               conv_tol_grad=getattr(mf, "conv_tol_grad", None),
               grad_norm=float(np.linalg.norm(g)) if g is not None else None)
    return mol, mf, scf


def make_td(mf, nstates):
    td = mf.TDDFT(equilibrium_solvation=False)
    td.singlet = True
    td.nstates = nstates
    return td


def solver_space(td, nstates):
    """Search-space parameters real_eig derives (pyscf.tdscf._lr_eig.real_eig, installed version)."""
    from pyscf.tdscf import _lr_eig
    nocc = int((td._scf.mo_occ > 0).sum())
    nvir = int(td._scf.mo_occ.size - nocc)
    A = nocc * nvir
    inc = nstates if _lr_eig.MAX_SPACE_INC is None else max(nstates, min(_lr_eig.MAX_SPACE_INC, A // 2))
    ms = int(td.max_memory * 1e6 / 8 / (4 * A) / 2 - inc)
    warn = ms < nstates * 4 < A
    if warn:
        ms = inc * 2
    ms = min(max(ms, nstates * 4), A)
    return dict(nocc=nocc, nvir=nvir, A_size=A, space_inc=inc, max_space=ms, memory_warning=bool(warn),
                holder_mb=4 * A * ms * 8 / 1e6, max_memory_mb=td.max_memory,
                MAX_SPACE_INC=_lr_eig.MAX_SPACE_INC)


# --------------------------------------------------------------------------- residuals and tracking
def tolerant_map(saved, mf, vecs):
    """Rotate saved MO-basis RPA vectors into the current orbitals WITHOUT the production orbital-energy check
    (diagnostic only). Returns vectors and the orbital diagnostics."""
    S = mf.mol.intor_symmetric("int1e_ovlp")
    C0, occ0, e0 = saved["mo_coeff_0"], saved["mo_occ_0"], saved["mo_energy_0"]
    C1, occ1, e1 = mf.mo_coeff, mf.mo_occ, mf.mo_energy
    O = C0.T @ S @ C1
    o, v = occ0 > 0, occ0 == 0
    Uoo, Uvv = O[np.ix_(o, o)], O[np.ix_(v, v)]
    diag = dict(max_dmo_energy=float(np.max(np.abs(e0 - e1))),
                occ_sv=[float(x) for x in np.linalg.svd(Uoo, compute_uv=False)[[0, -1]]],
                vir_sv=[float(x) for x in np.linalg.svd(Uvv, compute_uv=False)[[0, -1]]])
    no, nv = int(o.sum()), int(v.sum())
    n = no * nv
    out = []
    for vec in vecs:
        x = (Uoo.T @ vec[:n].reshape(no, nv) @ Uvv).ravel()
        y = (Uoo.T @ vec[n:].reshape(no, nv) @ Uvv).ravel()
        out.append(np.hstack([x, y]))
    return np.array(out), diag


def normalise(Z):
    n = Z.shape[1] // 2
    out = []
    for z in Z:
        m = z[:n] @ z[:n] - z[n:] @ z[n:]
        out.append(z / np.sqrt(abs(m)) * (1 if m > 0 else np.nan))
    return np.array(out)


def residuals(td, Z, omega):
    """True residual norms (real_eig convention) and Rayleigh-type quotient for each normalised vector."""
    vind, _ = td.gen_vind(td._scf)
    Z = normalise(Z)
    n = Z.shape[1] // 2
    ax = np.asarray(vind(Z))
    u1, u2 = ax[:, :n], -ax[:, n:]                      # A X + B Y,  B X + A Y
    X, Y = Z[:, :n], Z[:, n:]
    rx = u1 - omega[:, None] * X
    ry = u2 + omega[:, None] * Y
    r = np.sqrt((rx ** 2).sum(1) + (ry ** 2).sum(1))
    # quotient omega_q = (X.(AX+BY) + Y.(BX+AY)) / (X.X - Y.Y) for comparison with the stored Ritz value
    q = (X * u1).sum(1) + (Y * u2).sum(1)
    return r, q


def metric_overlap(Za, Zb):
    """|<a|b>| in the RPA metric (X.X' - Y.Y') for normalised vectors."""
    n = Za.shape[1] // 2
    return np.abs(Za[:, :n] @ Zb[:, :n].T - Za[:, n:] @ Zb[:, n:].T)


def clusters(e_au, gap_ev=0.05):
    e = np.asarray(e_au) * 27.211386
    out, cur = [], [0]
    for i in range(1, len(e)):
        if e[i] - e[i - 1] < gap_ev:
            cur.append(i)
        else:
            out.append(cur)
            cur = [i]
    out.append(cur)
    return [c for c in out if len(c) > 1]


def cmd_residuals(a):
    work = os.path.join(DIAG, a.tag)
    mol, mf, scf = build(a.system, a.func, work)
    td = make_td(mf, a.nstates)
    rows, prev = [], None
    report = dict(system=a.system, func=a.func, scf=scf, checkpoints=[])
    for path in a.ckpts:
        d = dict(np.load(path, allow_pickle=True))
        meta = json.loads(str(d["meta"]))
        Z, odiag = tolerant_map(d, mf, d["vectors"])
        Z = normalise(Z)
        om = np.asarray(d["energies"], float)
        t0 = time.time()
        r, q = residuals(td, Z, om)
        ent = dict(path=os.path.relpath(path, HERE), sha256=sha256(path), total_cycles=meta.get("total_cycles"),
                   converged_tol=meta.get("converged_tol"), updated=meta.get("updated"),
                   orbital_mapping=odiag, seconds=time.time() - t0,
                   roots=[dict(k=k, e_ev=float(om[k] * 27.211386), resid=float(r[k]),
                               quotient_minus_ritz_ev=float((q[k] - om[k]) * 27.211386),
                               saved_conv=bool(d["conv"][k])) for k in range(len(om))],
                   clusters_0p05ev=clusters(om))
        if prev is not None:
            S = metric_overlap(prev[0], Z)
            best = S.argmax(1)
            ent["tracking_vs_previous"] = [dict(prev_k=i, best_k=int(best[i]), overlap=float(S[i, best[i]]),
                                                second=float(np.sort(S[i])[-2])) for i in range(len(best))]
            # principal angles of each current cluster's subspace vs the previous same-index subspace
            angs = []
            for c in clusters(om):
                Sa = S[np.ix_(c, c)]
                sv = np.linalg.svd(Sa, compute_uv=False)
                angs.append(dict(cluster=c, min_cos=float(sv.min())))
            ent["cluster_subspace_vs_previous"] = angs
        report["checkpoints"].append(ent)
        prev = (Z, om)
    out = os.path.join(work, "residuals.json")
    json.dump(report, open(out, "w"), indent=1)
    print(json.dumps({os.path.basename(c["path"]): dict(cycles=c["total_cycles"],
                      unconverged_0p01=[r["k"] for r in c["roots"] if r["resid"] > 1e-2],
                      max_resid=max(r["resid"] for r in c["roots"])) for c in report["checkpoints"]}, indent=1))
    print(out)


# --------------------------------------------------------------------------- segments
ROOT_RE = re.compile(r"root (\d+)\s+(?:converged\s+)?\|r\|=\s*([0-9.eE+-]+)\s+e=")


class CountingVind:
    """Counts vectors passed through the response operator; each call is also appended to calls.jsonl at once,
    so a segment stopped early (bounded compute) still has its operator-application record."""
    stream = None

    def __init__(self, vind):
        self.vind, self.calls = vind, []

    def __call__(self, x):
        t0 = time.time()
        out = self.vind(x)
        rec = dict(nvec=int(np.asarray(x).shape[0]), seconds=time.time() - t0, t_end=time.time())
        self.calls.append(rec)
        if CountingVind.stream is not None:
            CountingVind.stream.write(json.dumps(rec) + "\n")
            CountingVind.stream.flush()
        return out


def cmd_segment(a):
    tag = a.tag
    work = os.path.join(DIAG, tag)
    if os.path.exists(os.path.join(work, "segment.json")):
        print(f"{work}/segment.json exists; refusing to overwrite")
        return
    os.makedirs(work, exist_ok=True)
    # private SCF directory per segment (concurrent segments must not share a chkfile); seeded from the shared
    # diagnostic SCF copy when one exists, otherwise from the production chk (build() copies it)
    scfdir = os.path.join(work, "scf")
    os.makedirs(scfdir, exist_ok=True)
    shared = os.path.join(DIAG, f"scf_{a.system}_{a.func}", f"scf_{a.func}.chk")
    if os.path.isfile(shared) and not os.path.isfile(os.path.join(scfdir, f"scf_{a.func}.chk")):
        shutil.copy2(shared, os.path.join(scfdir, f"scf_{a.func}.chk"))
    mol, mf, scf = build(a.system, a.func, scfdir, max_memory=a.max_memory)
    cdir = os.path.join(work, "td_ckpt")
    total0 = 0
    if a.ckpt:
        os.makedirs(cdir, exist_ok=True)
        shutil.copy2(a.ckpt, os.path.join(cdir, "latest.npz"))
        src_fp = os.path.join(os.path.dirname(a.ckpt), "fingerprint.json")
        if os.path.isfile(src_fp):
            shutil.copy2(src_fp, os.path.join(cdir, "fingerprint.json"))
        total0 = int(json.loads(str(np.load(a.ckpt, allow_pickle=True)["meta"]))["total_cycles"])
    td = make_td(mf, a.nstates)
    td.max_memory = a.max_memory or td.max_memory
    log_path = os.path.join(work, "solver.log")
    td.stdout = open(log_path, "w")
    CountingVind.stream = open(os.path.join(work, "calls.jsonl"), "w")
    td.verbose = 6
    counters = []
    orig = td.gen_vind

    def gen_vind(mf_=None):
        vind, hdiag = orig(mf_)
        c = CountingVind(vind)
        counters.append(c)
        return c, hdiag
    td.gen_vind = gen_vind
    lines = []
    t0 = time.time()
    r = TC.solve(td, mf, tda=False, nstates=a.nstates, conv_tol=a.conv_tol, ckpt_dir=cdir, singlet=True,
                 equilibrium_solvation=False, chunk=a.chunk, max_total_cycles=total0 + a.cycles,
                 log=lambda *m: lines.append(" ".join(map(str, m))))
    wall = time.time() - t0
    td.stdout.close()
    # per-cycle per-root residuals: each solver cycle logs one line per root
    resid, cur = [], {}
    for ln in open(log_path):
        m = ROOT_RE.search(ln)
        if m:
            k = int(m.group(1))
            if k in cur:
                resid.append([cur[i] for i in sorted(cur)])
                cur = {}
            cur[k] = float(m.group(2))
    if cur:
        resid.append([cur[i] for i in sorted(cur)])
    calls = [c for cc in counters for c in cc.calls]
    out = dict(system=a.system, func=a.func, nstates=a.nstates, chunk=a.chunk, cycles_budget=a.cycles,
               conv_tol=a.conv_tol, max_memory=td.max_memory, space=solver_space(td, a.nstates), scf=scf,
               start_ckpt=None if not a.ckpt else dict(path=os.path.relpath(a.ckpt, HERE), sha256=sha256(a.ckpt),
                                                       total_cycles=total0),
               status=r["status"], total_cycles=r["total_cycles"], wall_seconds=wall, ckpt_log=lines,
               operator_calls=calls, operator_vectors=int(sum(c["nvec"] for c in calls)),
               operator_seconds=float(sum(c["seconds"] for c in calls)),
               residual_per_cycle=resid, final_energies_ev=(np.asarray(r["energies_au"]) * 27.211386).tolist(),
               final_conv=np.asarray(r["conv"]).tolist())
    json.dump(out, open(os.path.join(work, "segment.json"), "w"), indent=1)
    last = resid[-1] if resid else []
    print(json.dumps(dict(tag=tag, status=r["status"], cycles=len(resid), wall_s=round(wall),
                          operator_vectors=out["operator_vectors"],
                          unconverged=[k for k, v in enumerate(last) if v > a.conv_tol_stage],
                          max_resid_last=max(last) if last else None), indent=1))


def cmd_stability(a):
    """Internal (real, singlet) stability of the SCF solution used by the response calculation. An internal
    instability (negative lowest eigenvalue of the orbital Hessian) would invalidate RPA and could cause
    erratic response convergence."""
    from pyscf.scf import stability
    work = os.path.join(DIAG, f"stab_{a.system}_{a.func}")
    mol, mf, scf = build(a.system, a.func, work)
    t0 = time.time()
    lines = []
    mf.verbose, mf.stdout = 4, open(os.path.join(work, "stability.log"), "w")
    mo_i, mo_e, stable_i, stable_e = mf.stability(internal=True, external=False, return_status=True)
    mf.stdout.close()
    for ln in open(os.path.join(work, "stability.log")):
        if "stability" in ln.lower() or "lowest eigenvalue" in ln.lower():
            lines.append(ln.strip())
    out = dict(system=a.system, func=a.func, scf=scf, internal_stable=bool(stable_i), seconds=time.time() - t0,
               log=lines)
    json.dump(out, open(os.path.join(work, "stability.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


def cmd_config(a):
    work = os.path.join(DIAG, f"scf_{a.system}_{a.func}")
    mol, mf, scf = build(a.system, a.func, work)
    td = make_td(mf, a.nstates)
    import pyscf
    from pyscf.tdscf import _lr_eig
    commit = subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    cfg = dict(commit=commit, pyscf=pyscf.__version__, numpy=np.__version__, python=sys.version.split()[0],
               system=a.system, functional=a.func, xc=mf.xc, basis=mol.basis, cartesian=mol.cart,
               natm=mol.natm, nelectron=mol.nelectron, nao=mol.nao,
               density_fit="def2-universal-jkfit", grid_level=mf.grids.level,
               solvent=dict(method=mf.with_solvent.method, eps=mf.with_solvent.eps, lebedev=mf.with_solvent.lebedev_order,
                            td_equilibrium_solvation=False),
               td_class=type(td).__mro__[0].__name__, td_mro=[c.__name__ for c in type(td).__mro__[:6]],
               response="full TD-DFT (RPA, non-Hermitian real_eig), singlets", nstates=a.nstates,
               solver=f"{_lr_eig.real_eig.__module__}.real_eig", solver_file=_lr_eig.__file__,
               solver_sha256=sha256(_lr_eig.__file__),
               convergence="per-root residual norm <= conv_tol (vectors normalised X.X-Y.Y=1); staged 1e-2 -> 1e-5",
               td_lindep=td.lindep, td_level_shift=td.level_shift, td_deg_eia_thresh=td.deg_eia_thresh,
               space=solver_space(td, a.nstates), scf=scf,
               production=dict(chunk=2, max_total_cycles=400, threads=2, mem_reservation_mb=3000,
                               restart_start_set="30 Ritz vectors + 30 standard initial-guess vectors, SVD-orthonormalised"))
    out = os.path.join(DIAG, f"config_{a.system}_{a.func}.json")
    json.dump(cfg, open(out, "w"), indent=1, default=str)
    print(json.dumps(cfg, indent=1, default=str))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("config")
    c.add_argument("--system", default="L1")
    c.add_argument("--func", default="CAM-B3LYP")
    c.add_argument("--nstates", type=int, default=30)
    r = sub.add_parser("residuals")
    r.add_argument("--system", default="L1")
    r.add_argument("--func", required=True)
    r.add_argument("--nstates", type=int, default=30)
    r.add_argument("--tag", required=True)
    r.add_argument("ckpts", nargs="+")
    s = sub.add_parser("segment")
    s.add_argument("--system", default="L1")
    s.add_argument("--func", required=True)
    s.add_argument("--nstates", type=int, default=30)
    s.add_argument("--ckpt", default=None, help="checkpoint to COPY as the start (omit: fresh start)")
    s.add_argument("--chunk", type=int, required=True)
    s.add_argument("--cycles", type=int, required=True, help="Davidson-cycle budget of this segment")
    s.add_argument("--conv-tol", type=float, default=1e-5, help="final target (stages 1e-2 ... target)")
    s.add_argument("--conv-tol-stage", type=float, default=1e-2, help="for the printed summary only")
    s.add_argument("--max-memory", type=int, default=None, help="td.max_memory (sets the real_eig max_space)")
    s.add_argument("--tag", required=True)
    st = sub.add_parser("stability")
    st.add_argument("--system", default="L1")
    st.add_argument("--func", required=True)
    a = p.parse_args(argv)
    qc.set_threads()
    dict(config=cmd_config, residuals=cmd_residuals, segment=cmd_segment, stability=cmd_stability)[a.cmd](a)


if __name__ == "__main__":
    main()
