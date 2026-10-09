"""Fe(III) increment analysis: validated inputs, shared stick-based broadening, neutral naming."""

import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import completion  # noqa: E402
import fe_increment as FI  # noqa: E402


def _run(d, level, e, f, basis="6-31g*", marker=True, geom="converged"):
    d.mkdir(parents=True)
    summ = dict(level=level, basis=basis, solvent="pcm", eps=78.36, tda=True, cartesian_d=True,
                tddft={"B3LYP": {}}, geometry_status=dict(status=geom))
    (d / "summary.json").write_text(json.dumps(summ))
    st = d / f"{level}_B3LYP_states.csv"
    pd.DataFrame(dict(State=np.arange(1, len(e) + 1), Wavelength_nm=1239.84 / np.array(e), Energy_eV=e,
                      Oscillator_Strength=f)).to_csv(st, index=False)
    if marker:
        completion.write_marker(str(d), summ, [str(st), str(d / "summary.json")], True, dict(status=geom))
    return str(d)


def test_stick_rebuild_settings_check_and_names(tmp_path, monkeypatch):
    lig = _run(tmp_path / "lig", "level1", [3.2, 4.4, 5.0, 6.5], [0.1, 0.4, 0.3, 0.5])
    cpx = _run(tmp_path / "cpx", "level3_catecholate", [2.0, 3.1, 4.4, 6.5], [0.02, 0.1, 0.4, 0.5])
    eps, st, level, meta = FI.eps_of(lig, "B3LYP")
    assert np.all(np.isfinite(eps)) and FI.GRID[0] == 280.0       # covers 280 nm from the sticks
    with pytest.raises(ValueError):                                # basis differs from the ligand
        FI.eps_of(_run(tmp_path / "x", "level3_tropolonate", [2.0], [0.1], basis="def2-svp"), "B3LYP",
                  reference=meta["settings"])
    with pytest.raises(ValueError):                                # no completion marker
        FI.eps_of(_run(tmp_path / "y", "level3_tropolonate", [2.0], [0.1], marker=False), "B3LYP")
    with pytest.raises(ValueError):                                # provisional not accepted by default
        FI.eps_of(_run(tmp_path / "z", "level3_tropolonate", [2.0], [0.1], geom="exploratory_xtb"), "B3LYP")
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["x"])
    FI.main(["--ligand", lig, "--complex", cpx, str(tmp_path / "z"), "--outdir", str(out), "--allow-provisional"])
    tab = pd.read_csv(out / "fe_increment_metrics.csv")
    assert "strongest_vis_state_nm" in tab and not any("lmct" in c for c in tab.columns)
    assert np.isclose(tab.strongest_vis_state_nm.iloc[1], 1239.84 / 2.0, rtol=1e-3)
    assert bool(tab.exploratory.iloc[2]) and tab.geometry.iloc[2] == "exploratory_xtb"
