"""Completion markers and run_phase1 end-to-end status handling (tiny molecule, real PySCF runs)."""

import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import completion  # noqa: E402
import molecules  # noqa: E402
import run_phase1  # noqa: E402

XYZ = "3\n\nO 0 0 0\nH 0 0.76 0.59\nH 0 -0.76 0.59\n"


@pytest.fixture
def tiny(monkeypatch):
    monkeypatch.setitem(molecules.MOLECULES, "tiny", dict(name="water", formula="H2O", molar_mass=18.015,
                                                          charge=0, multiplicity=1))
    return "tiny"


def _run(tmp_path, level, extra=()):
    xyz = tmp_path / "w.xyz"
    xyz.write_text(XYZ)
    out = tmp_path / "out"
    args = ["--level", level, "--basis", "6-31g", "--skip-opt", "--start-xyz", str(xyz), "--solvent", "none",
            "--nstates", "3", "--outdir", str(out), "--verbose", "0", "--max-memory", "1000", *extra]
    return run_phase1.main(args), out, xyz


def test_td_problems_detects_each_failure():
    ok = completion.td_problems([1, 2, 3], [0.1, 0, 0.2], [True] * 3, 3)
    assert ok == []
    assert completion.td_problems([1, 2], [0.1, 0.2], [True, True], 3)
    assert completion.td_problems([1, 2, 3], [0.1, 0.2, 0.3], [True, False, True], 3)
    assert completion.td_problems([1, np.nan, 3], [0.1, 0.2, 0.3], [True] * 3, 3)
    assert completion.td_problems([1, 2, 3], [0.1, -0.2, 0.3], [True] * 3, 3)
    assert completion.td_problems([-1, 2, 3], [0.1, 0.2, 0.3], [True] * 3, 3)


def test_unverified_start_geometry_is_provisional_and_tamper_detected(tiny, tmp_path):
    s, out, _ = _run(tmp_path, tiny)
    assert s["completion_status"] == "provisional_geometry" and run_phase1._exit_code(s) == 3
    ok, st, why = completion.validate_run(str(out))
    assert not ok and st == "provisional_geometry"
    ok, _, _ = completion.validate_run(str(out), allow_provisional=True)
    assert ok
    with open(out / "tiny_B3LYP_spectrum.csv", "a") as fh:            # modify an artifact afterwards
        fh.write("\n")
    ok, _, why = completion.validate_run(str(out), allow_provisional=True)
    assert not ok and any("changed" in w for w in why)


def test_accepted_geometry_sidecar_gives_production(tiny, tmp_path):
    xyz = tmp_path / "w.xyz"
    xyz.write_text(XYZ)
    (tmp_path / "w.xyz.geometry.json").write_text(json.dumps(dict(status="converged")))
    s, out, _ = _run(tmp_path, tiny)
    assert s["completion_status"] == "production" and run_phase1._exit_code(s) == 0
    rcs = s["tddft"]["B3LYP"]["root_count_sensitivity"]
    assert [r["window_nm"] for r in rcs] == [[260.0, 750.0]] * 2 and [r["fwhm_ev"] for r in rcs] == [0.3, 0.6]
    ok, st, _ = completion.validate_run(str(out), require=dict(nstates=3, tda=False))
    assert ok and st == "production"
    ok, _, why = completion.validate_run(str(out), require=dict(nstates=30))
    assert not ok and any("nstates" in w for w in why)


def test_unconverged_roots_are_diagnostic_only(tiny, tmp_path):
    s, out, _ = _run(tmp_path, tiny, ["--td-conv-tol", "1e-12", "--td-chunk", "1", "--td-max-cycles", "1"])
    assert s["completion_status"] == "diagnostic_only" and run_phase1._exit_code(s) == 4
    assert not os.path.exists(out / completion.MARKER) and os.path.exists(out / completion.DIAG)
    ok, st, _ = completion.validate_run(str(out), allow_provisional=True)
    assert not ok and st == "diagnostic_only"


def test_partial_output_without_marker_is_rejected(tmp_path):
    (tmp_path / "summary.json").write_text("{}")                       # nonempty summary alone
    ok, st, _ = completion.validate_run(str(tmp_path), allow_provisional=True)
    assert not ok and st == "missing"
