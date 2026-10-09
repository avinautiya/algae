"""Regression tests: XYZ parsing (blank comments, frames, malformed/truncated input) and the geometry
optimiser's history/finality handling."""

import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import qc  # noqa: E402

W1 = "3\n\nO 0 0 0\nH 0 0.76 0.59\nH 0 -0.76 0.59\n"          # blank comment line is valid XYZ
W2 = "3\nframe 2\nO 0 0 0.01\nH 0 0.75 0.59\nH 0 -0.75 0.59\n"


def _write(tmp_path, text, name="g.xyz"):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_blank_comment_single_frame(tmp_path):
    frames, w = qc.read_xyz_frames(_write(tmp_path, W1))
    assert len(frames) == 1 and frames[0][0] == "" and len(frames[0][1]) == 3 and not w


def test_multiple_frames_last_returned(tmp_path):
    p = _write(tmp_path, W1 + W2)
    frames, _ = qc.read_xyz_frames(p)
    assert len(frames) == 2
    assert qc.xyz_file_to_atom_block(p).splitlines()[0].split()[3] == "0.0100000000"


def test_malformed_coordinates_raise(tmp_path):
    with pytest.raises(qc.XYZError):
        qc.read_xyz_frames(_write(tmp_path, "3\nc\nO 0 0 x\nH 0 1 0\nH 0 -1 0\n"))


def test_truncated_last_frame_recovers_previous(tmp_path):
    frames, w = qc.read_xyz_frames(_write(tmp_path, W1 + W2 + "3\ncut\nO 0 0 0\n"))
    assert len(frames) == 2 and w and "truncated" in w[0]


def test_truncated_only_frame_raises(tmp_path):
    with pytest.raises(qc.XYZError):
        qc.read_xyz_frames(_write(tmp_path, "3\nc\nO 0 0 0\n"))


def test_malformed_after_valid_frame_recovers(tmp_path):
    frames, w = qc.read_xyz_frames(_write(tmp_path, W1 + "3\nc\nO 0 0 0\nH a b c\nH 0 1 0\n"))
    assert len(frames) == 1 and w


def test_no_recovery_when_strict(tmp_path):
    with pytest.raises(qc.XYZError):
        qc.read_xyz_frames(_write(tmp_path, W1 + "3\ncut\nO 0 0 0\n"), recover=False)


def test_optimizer_final_only_when_converged_and_history_appended(tmp_path):
    from pyscf import gto
    mol = gto.M(atom="O 0 0 0; H 0 0.80 0.62; H 0 -0.80 0.62", basis="sto-3g", verbose=0)
    _, _, _, conv1 = qc.optimize_geometry(mol, "B3LYP", None, maxsteps=1, workdir=str(tmp_path), tag="t")
    assert not conv1
    assert not os.path.exists(tmp_path / "t_final.xyz")
    assert os.path.exists(tmp_path / "t_unconverged_endpoint.xyz")
    rec = json.load(open(tmp_path / "t_opt_record.json"))
    assert rec["converged"] is False and set(rec["thresholds"]) == set(qc.GAUSSIAN_CRITERIA)
    n1 = len(open(tmp_path / "t_steps.jsonl").read().splitlines())
    mol2 = gto.M(atom=qc.xyz_file_to_atom_block(str(tmp_path / "t_last.xyz")), basis="sto-3g", verbose=0)
    _, _, _, conv2 = qc.optimize_geometry(mol2, "B3LYP", None, maxsteps=50, workdir=str(tmp_path), tag="t")
    assert conv2 and os.path.exists(tmp_path / "t_final.xyz")
    rows = [json.loads(l) for l in open(tmp_path / "t_steps.jsonl")]
    assert len(rows) > n1 and {r["segment"] for r in rows} == {1, 2}       # history kept, segments recorded
    frames, _ = qc.read_xyz_frames(str(tmp_path / "t_traj.xyz"))
    assert len(frames) == len(rows)
