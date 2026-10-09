"""Cached results are reused only on matching content fingerprints; MCMC chains are reproducible."""

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
for d in ("..", "../../phase3", "../../phase4", "../../phase1"):
    sys.path.insert(0, os.path.join(HERE, d))

import provenance as PV  # noqa: E402


def test_fingerprint_sensitive_to_content_not_identity():
    a = dict(x=np.arange(5.0), s="a", d={"k": 1})
    b = dict(d={"k": 1}, s="a", x=np.arange(5.0))
    assert PV.fingerprint(a) == PV.fingerprint(b)
    c = dict(a, x=np.arange(5.0) + 1e-12)
    assert PV.fingerprint(a) != PV.fingerprint(c)
    assert PV.fingerprint(np.arange(3, dtype=np.float32)) != PV.fingerprint(np.arange(3, dtype=np.float64))


def test_run_or_load_rejects_same_row_count_with_other_inputs(tmp_path, monkeypatch):
    import forward_model
    import run_phase3
    calls = []

    def fake_eval(kw, design, workers=1):
        calls.append(len(design))
        return pd.DataFrame(dict(y=[d["p"] * kw["scale"] for d in design]))
    monkeypatch.setattr(forward_model, "evaluate_many", fake_eval)
    path = str(tmp_path / "runs.csv")
    d1 = [dict(p=1.0), dict(p=2.0)]
    d2 = [dict(p=3.0), dict(p=4.0)]                           # same number of rows, different design
    kw = dict(scale=2.0, phase1_l2=None, biosnicar=None)
    r1 = run_phase3.run_or_load(path, d1, kw, 1, True, "t")
    r1b = run_phase3.run_or_load(path, d1, kw, 1, True, "t")
    assert calls == [2] and r1b.y.tolist() == r1.y.tolist()        # reused
    r2 = run_phase3.run_or_load(path, d2, kw, 1, True, "t")
    assert calls == [2, 2] and r2.y.tolist() == [6.0, 8.0]         # old code reused r1 here
    run_phase3.run_or_load(path, d2, dict(kw, scale=3.0), 1, True, "t")
    assert len(calls) == 3                                        # settings change -> recompute
    with open(path, "a") as fh:                                   # result edited after writing
        fh.write("1\n")
    run_phase3.run_or_load(path, d2, dict(kw, scale=3.0), 1, True, "t")
    assert len(calls) == 4


def test_cached_calibration_keys_on_settings_and_content(monkeypatch):
    import cell_optics as co
    import tddft_calibration as TC
    n = []
    monkeypatch.setattr(TC, "calibrate", lambda spec, **kw: n.append(kw) or object())
    TC._CACHE.clear()
    s = co.MolecularSpectrum("x", 300.0, 0.3, np.array([2.5, 3.0]), np.array([0.1, 0.2]), source="same")
    TC.cached_calibration(s, seed=0)
    TC.cached_calibration(s, seed=0)
    TC.cached_calibration(s, seed=1)                                # old key ignored kw
    TC.cached_calibration(co.MolecularSpectrum("x", 310.0, 0.3, s.energies_ev, s.osc, source="same"), seed=0)
    assert len(n) == 3


def test_phase1_fingerprint_follows_content(tmp_path):
    (tmp_path / "level2_B3LYP_states.csv").write_text("a\n1\n")
    (tmp_path / "summary.json").write_text("{}")
    f1 = PV.phase1_fingerprint(str(tmp_path))
    (tmp_path / "level2_B3LYP_states.csv").write_text("a\n2\n")
    assert PV.phase1_fingerprint(str(tmp_path)) != f1


def test_calibration_chains_reproducible():
    import cell_optics as co
    import tddft_calibration as TC
    s = co.MolecularSpectrum("x", 474.4, 0.3, np.array([2.95, 3.25, 3.72, 4.19, 4.6, 5.1]),
                             np.array([0.07, 0.15, 0.2, 0.27, 0.1, 0.3]), source="test")
    a = TC.calibrate(s, n_walkers=12, n_steps=120, burn=40, n_stage2=6, draws_per=3, seed=3, verbose=False)
    b = TC.calibrate(s, n_walkers=12, n_steps=120, burn=40, n_stage2=6, draws_per=3, seed=3, verbose=False)
    np.random.seed(12345)                                         # global state must not matter
    c = TC.calibrate(s, n_walkers=12, n_steps=120, burn=40, n_stage2=6, draws_per=3, seed=3, verbose=False)
    assert np.array_equal(a.samples, b.samples) and np.array_equal(a.samples, c.samples)


def test_emulator_save_is_atomic(tmp_path):
    import emulator as E
    em = E.Emulator({"log_b": np.array([1.0, 2.0])}, {"bands": np.zeros((2, 4))}, dict(tag="t"))
    p = str(tmp_path / "e.npz")
    em.save(p)
    assert os.listdir(tmp_path) == ["e.npz"] and E.Emulator.load(p).meta["tag"] == "t"


def test_fit_magnitude_has_no_silent_covariance_fallback(monkeypatch):
    import pytest
    import cell_optics as co
    import tddft_calibration as TC
    from scipy import optimize
    s = co.MolecularSpectrum("x", 474.4, 0.3, np.array([2.95, 3.25, 3.72, 4.19]), np.array([0.07, 0.15, 0.2, 0.27]))
    wl_h, S, sS, wl_m, E, sE, D = TC._data()
    real = optimize.minimize

    def failing(*a, **k):
        r = real(*a, **k)
        r.success = False
        return r
    monkeypatch.setattr(optimize, "minimize", failing)
    with pytest.raises(RuntimeError):
        TC._fit_magnitude(s, 0.0, 0.5, wl_m, E, sE, D)
    _, ok = TC._fit_magnitude(s, 0.0, 0.5, wl_m, E, sE, D, return_status=True)
    assert ok is False


def test_stage2_no_boundary_pileup_and_quadrature_refinement():
    import cell_optics as co
    import tddft_calibration as TC
    s = co.MolecularSpectrum("x", 474.4, 0.3, np.array([2.95, 3.25, 3.72, 4.19, 4.6, 5.1]),
                             np.array([0.07, 0.15, 0.2, 0.27, 0.1, 0.3]))
    wl_h, S, sS, wl_m, E, sE, D = TC._data()
    rng = np.random.default_rng(0)
    out = {}
    for n in (1, 2):
        z, dg = TC._stage2_grid(s, 0.05, 0.6, wl_m, E, sE, D, rng, 40000, n_lf=81 * n, n_phi=101 * n - (n - 1),
                                n_ls=56 * n, stride=2, rho=(0.9, 0.99))
        assert np.all((z[:, 1] >= 0) & (z[:, 1] <= 1))
        assert np.mean(z[:, 1] == 0.0) == 0 and np.mean(z[:, 1] == 1.0) == 0      # no clipping pile-up
        out[n] = z.mean(0), z.std(0), dg["rho_draws"]
    # ln f and phi (which enter the optics) converge to < 0.01; the discrepancy ln s (not used downstream)
    # to < 0.05 at the default resolution (records/repair_stage2_refinement.json has the full study)
    assert np.allclose(out[1][0][:2], out[2][0][:2], atol=0.01) and abs(out[1][0][2] - out[2][0][2]) < 0.05
    assert np.allclose(out[1][1][:2], out[2][1][:2], rtol=0.1, atol=0.005)
    assert set(np.unique(out[1][2])) <= {0.9, 0.99}
