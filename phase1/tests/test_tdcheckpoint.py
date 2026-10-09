"""Regression tests for restartable TD/TDA solves (phase1/tdcheckpoint.py).

These verify restart MECHANICS on a small solvated molecule (formaldehyde, PCM). They do not establish
the accuracy of any approximation for the target pigment.
"""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import qc  # noqa: E402
import tdcheckpoint as TC  # noqa: E402

GEOM = "C 0 0 0; O 0 0 1.21; H 0 0.94 -0.59; H 0 -0.94 -0.59"
GEOM2 = "C 0 0 0; O 0 0 1.23; H 0 0.94 -0.59; H 0 -0.94 -0.59"
NST = 4


def _mf(geom=GEOM, xc="B3LYP", basis="6-31g*", solvent="pcm"):
    from pyscf import gto
    mol = gto.M(atom=geom, basis=basis, cart=True, verbose=0)
    mf = qc.make_mf(mol, xc, solvent)
    mf.kernel()
    assert mf.converged
    return mf


def _td(mf, tda=True):
    td = (mf.TDA(equilibrium_solvation=False) if tda else mf.TDDFT(equilibrium_solvation=False)) \
        if getattr(mf, "with_solvent", None) is not None else (mf.TDA() if tda else mf.TDDFT())
    td.singlet, td.nstates = True, NST
    return td


def _solve(mf, d, tol=1e-5, tda=True, **kw):
    td = _td(mf, tda)
    r = TC.solve(td, mf, tda=tda, nstates=NST, conv_tol=tol, ckpt_dir=str(d), singlet=True, log=lambda *a: None, **kw)
    f = np.asarray(td.oscillator_strength(gauge="length")) if r["status"] in ("converged",) else None
    return r, td, f


@pytest.fixture(scope="module")
def mf():
    return _mf()


@pytest.fixture(scope="module")
def reference(mf):
    td = _td(mf)
    td.conv_tol = 1e-5
    td.kernel()
    return np.asarray(td.e), np.asarray(td.oscillator_strength(gauge="length"))


@pytest.mark.parametrize("tda", [True, False])
def test_uninterrupted_matches_plain(mf, tmp_path, tda):
    td = _td(mf, tda)
    td.conv_tol = 1e-5
    td.kernel()
    r, td2, f = _solve(mf, tmp_path, tda=tda)
    assert r["status"] == "converged"
    assert np.allclose(td2.e, td.e, atol=1e-7)
    assert np.allclose(f, td.oscillator_strength(gauge="length"), atol=1e-5)


def test_interrupted_before_stage_completion_then_resumed(mf, tmp_path, reference):
    r1, _, _ = _solve(mf, tmp_path, chunk=1, max_total_cycles=2)          # budget runs out mid-stage
    assert r1["status"] == "max_cycles_unconverged"
    r2, td, f = _solve(mf, tmp_path, chunk=3)
    assert r2["status"] == "converged" and r2["resumed"]
    assert np.allclose(td.e, reference[0], atol=1e-7)
    assert np.allclose(f, reference[1], atol=1e-5)


def test_between_stages_and_changed_final_tolerance(mf, tmp_path, reference):
    r1, _, _ = _solve(mf, tmp_path, tol=1e-3)
    assert r1["status"] == "converged" and r1["stage_history"][-1]["tol"] == 1e-3
    r2, td, f = _solve(mf, tmp_path, tol=1e-5)                              # tighter target reuses vectors
    assert r2["resumed"] and [h["tol"] for h in r2["stage_history"]][-1] == 1e-5
    assert np.allclose(td.e, reference[0], atol=1e-7)
    assert np.allclose(f, reference[1], atol=1e-5)


def test_max_cycle_without_convergence_is_not_success(mf, tmp_path):
    r, _, _ = _solve(mf, tmp_path, tol=1e-9, chunk=1, max_total_cycles=1)
    assert r["status"] == "max_cycles_unconverged"
    assert not r["conv"].all()


@pytest.mark.parametrize("change", ["geometry", "functional", "basis", "solvent"])
def test_changed_operator_same_nstates_is_not_reused(mf, tmp_path, change):
    _solve(mf, tmp_path, tol=1e-3)
    other = {"geometry": dict(geom=GEOM2), "functional": dict(xc="CAM-B3LYP"),
             "basis": dict(basis="6-31g"), "solvent": dict(solvent=None)}[change]
    r, _, _ = _solve(_mf(**other), tmp_path, tol=1e-3)
    assert not r["resumed"]


def test_changed_mo_phase_is_transformed(mf, tmp_path, reference):
    _solve(mf, tmp_path, tol=1e-3)
    mf2 = _mf()
    C = mf2.mo_coeff.copy()
    C[:, ::2] *= -1.0                                                       # arbitrary orbital sign flips
    mf2.mo_coeff = C
    r, td, f = _solve(mf2, tmp_path, tol=1e-5)
    assert r["resumed"] and r["status"] == "converged"
    assert np.allclose(td.e, reference[0], atol=1e-7)
    assert np.allclose(f, reference[1], atol=1e-5)


def test_noncanonical_rotation_invalidates(mf, tmp_path):
    """Rotating non-degenerate occupied orbitals changes PySCF's response operator (it uses orbital-energy
    differences), so such orbitals must not be accepted even though the occupied subspace is unchanged."""
    _solve(mf, tmp_path, tol=1e-3)
    saved = dict(np.load(os.path.join(tmp_path, "latest.npz")))
    mf2 = _mf()
    C = mf2.mo_coeff.copy()
    o = np.flatnonzero(mf2.mo_occ > 0)[-2:]
    c, s = np.cos(0.3), np.sin(0.3)
    C[:, o] = C[:, o] @ np.array([[c, -s], [s, c]])
    mf2.mo_coeff = C
    with pytest.raises(TC.CheckpointError):
        TC.map_vectors(dict(nspin=saved["nspin"], **{k: v for k, v in saved.items() if k.startswith("mo_")}),
                       mf2, saved["vectors"], True)


def test_changed_mo_subspace_invalidates(mf, tmp_path):
    _solve(mf, tmp_path, tol=1e-3)
    saved = dict(np.load(os.path.join(tmp_path, "latest.npz")))
    mf2 = _mf()
    C = mf2.mo_coeff.copy()
    i, a = np.flatnonzero(mf2.mo_occ > 0)[-1], np.flatnonzero(mf2.mo_occ == 0)[0]
    c, s = np.cos(0.2), np.sin(0.2)                                         # occupied-virtual mixing
    C[:, [i, a]] = C[:, [i, a]] @ np.array([[c, -s], [s, c]])
    mf2.mo_coeff = C
    with pytest.raises(TC.CheckpointError):
        TC.map_vectors(dict(nspin=saved["nspin"], **{k: v for k, v in saved.items() if k.startswith("mo_")}),
                       mf2, saved["vectors"], True)


def test_corrupt_latest_recovers_previous(mf, tmp_path):
    _solve(mf, tmp_path, tol=1e-3, chunk=1)                                 # several saves -> prev exists
    assert os.path.isfile(os.path.join(tmp_path, "latest.prev.npz"))
    with open(os.path.join(tmp_path, "latest.npz"), "wb") as fh:
        fh.write(b"garbage")
    r, _, _ = _solve(mf, tmp_path, tol=1e-3)
    assert r["resumed"]


def test_truncated_and_corrupt_pair_restarts_cleanly(mf, tmp_path, reference):
    _solve(mf, tmp_path, tol=1e-3, chunk=1)
    for name in ("latest.npz", "latest.prev.npz"):
        p = os.path.join(tmp_path, name)
        data = open(p, "rb").read()
        open(p, "wb").write(data[: len(data) // 3])                          # partial writes
    r, td, _ = _solve(mf, tmp_path)
    assert not r["resumed"] and r["status"] == "converged"
    assert any(n.endswith(tuple(f".corrupt.{x}" for x in [""])) or ".corrupt." in n for n in os.listdir(tmp_path))
    assert np.allclose(td.e, reference[0], atol=1e-7)


def test_already_completed_restart_verifies(mf, tmp_path, reference):
    _solve(mf, tmp_path)
    r, td, f = _solve(mf, tmp_path)
    assert r["resumed"] and r["status"] == "converged"
    assert np.allclose(td.e, reference[0], atol=1e-7)


def test_concurrent_writer_is_refused(mf, tmp_path):
    with TC.CheckpointStore(str(tmp_path)):
        with pytest.raises(TC.CheckpointError):
            with TC.CheckpointStore(str(tmp_path)):
                pass


def test_rounding_only_geometry_difference_reuses_vectors_without_stage_credit(mf, tmp_path, reference):
    """A geometry re-read from an .xyz rounds differently from the in-memory optimised one; a difference
    within GEOM_TOL_BOHR must not discard the checkpoint, but converged-stage credit is not inherited."""
    r1, _, _ = _solve(mf, tmp_path, tol=1e-3)
    assert r1["status"] == "converged"
    shifted = "C 0 0 0; O 0 0 1.210004; H 0 0.94 -0.59; H 0 -0.94 -0.59"        # 4e-6 A = 7.6e-6 Bohr
    r2, td, f = _solve(_mf(geom=shifted), tmp_path, tol=1e-5)
    assert r2["resumed"] and r2["rounding_match_bohr"] is not None
    assert 1e-6 < r2["rounding_match_bohr"] <= TC.GEOM_TOL_BOHR
    assert r2["status"] == "converged"
    assert np.allclose(r2["energies_au"], reference[0], atol=1e-4)     # same physics to the TD tolerance


def test_geometry_difference_beyond_tolerance_is_not_reused(mf, tmp_path):
    _solve(mf, tmp_path, tol=1e-3)
    shifted = "C 0 0 0; O 0 0 1.2101; H 0 0.94 -0.59; H 0 -0.94 -0.59"          # 1e-4 A = 1.9e-4 Bohr
    r, _, _ = _solve(_mf(geom=shifted), tmp_path, tol=1e-3)
    assert not r["resumed"]


def test_equivalent_except_rounding_unit():
    a = dict(xc="B3LYP", coords_bohr=[[0, 0, 0], [0, 0, 2.0]], nstates=4)
    assert TC.equivalent_except_rounding(a, dict(a, coords_bohr=[[0, 0, 0], [0, 0, 2.000001]]))[0]
    assert not TC.equivalent_except_rounding(a, dict(a, coords_bohr=[[0, 0, 0], [0, 0, 2.001]]))[0]
    assert not TC.equivalent_except_rounding(a, dict(a, xc="CAM-B3LYP"))[0]
    assert not TC.equivalent_except_rounding(a, dict(a, nstates=5))[0]
    assert not TC.equivalent_except_rounding(None, a)[0]
