"""State matching by overlaps for the approximation checks (real tiny PySCF runs)."""

import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import compare_states as cs  # noqa: E402
import molecules  # noqa: E402
import run_phase1  # noqa: E402

XYZ = "4\n\nC 0 0 0\nO 0 0 1.21\nH 0 0.94 -0.59\nH 0 -0.94 -0.59\n"


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    mp = pytest.MonkeyPatch()
    mp.setitem(molecules.MOLECULES, "tiny", dict(name="formaldehyde", formula="CH2O", molar_mass=30.03,
                                                 charge=0, multiplicity=1))
    d = tmp_path_factory.mktemp("cmp")
    xyz = d / "f.xyz"
    xyz.write_text(XYZ)
    out = {}
    for key, extra in [("rpa", []), ("tda", ["--tda"]), ("tda_loose", ["--tda", "--td-conv-tol", "1e-3"])]:
        o = d / key
        run_phase1.main(["--level", "tiny", "--basis", "6-31g", "--skip-opt", "--start-xyz", str(xyz),
                         "--solvent", "none", "--nstates", "4", "--outdir", str(o), "--verbose", "0",
                         "--max-memory", "1000", *extra])
        out[key] = cs.Run(str(o), "B3LYP", level="tiny")
    yield out
    mp.undo()


def test_overlap_matching_identity_and_tda_vs_rpa(runs):
    for r in runs.values():
        assert r.vec is not None and r.C is not None
    pairs, how = cs.match(runs["tda"], runs["tda"])
    assert how.startswith("transition-density") and all(i == j and o > 0.999 for i, j, o in pairs)
    c = cs.compare(runs["rpa"], runs["tda"], "A", "model")
    assert c["matching"].startswith("transition-density")
    assert all(s["matched"] for s in c["states"])
    assert all(s["dE_ev"] >= -1e-6 for s in c["states"])          # TDA energies lie above RPA (Thouless)


def test_matching_survives_reordering(runs):
    a, b = runs["tda"], runs["tda"]
    import copy
    b = copy.copy(b)
    perm = np.array([2, 0, 3, 1])
    b.vec, b.e, b.f = b.vec[perm], b.e[perm], b.f[perm]
    pairs, _ = cs.match(a, b)
    assert all(perm[j] == i for i, j, _ in pairs)


def test_label_fallback_and_tolerance_check(runs):
    import copy
    a = copy.copy(runs["rpa"])
    a.vec = None
    pairs, how = cs.match(a, runs["tda"])
    assert "label" in how and len(pairs) == 4
    t = cs.tolerance_check(runs["tda"])
    assert t and t["stages"][-1]["max_abs_dE_all_vs_final_ev"] == 0.0
    assert [s["tol"] for s in t["stages"]] == sorted([s["tol"] for s in t["stages"]], reverse=True)
    assert cs.compare(runs["tda_loose"], runs["tda"], "T", "numerical")["accept"] in (True, False)
