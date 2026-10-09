"""Held-out evaluation machinery: count likelihood with zeros, proper scoring, training-only settings."""

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))

import heldout as HO  # noqa: E402


def _df():
    return pd.DataFrame(dict(sample=["a", "b", "c"], cells=[1000.0, 0.0, 500.0], cells_counted=[16.0, 0.0, np.nan],
                             dataset=["s6_2017", "s6_2017", "sgris_2021"]))


def test_observation_likelihood_keeps_zeros_and_uses_counted_volume():
    L = HO.observation_loglik(_df())
    x = HO.XGRID
    # positive count: Poisson in N with V = N / cells -> maximum at B = cells
    assert abs(x[np.argmax(L[0])] - 3.0) < 0.011
    # zero count: P(N=0) = exp(-B V0), decreasing in B, ~1 for small B
    assert np.all(np.diff(L[1]) <= 1e-12) and L[1][0] > -0.01
    assert np.isclose(L[1][np.argmin(abs(x - 2))], -100 * HO.V_ZERO_ML)
    # log-normal observation error where the counted number is unknown
    assert abs(x[np.argmax(L[2])] - np.log10(500)) < 0.011


def test_log_score_is_normalised_density_score():
    L = np.zeros((1, HO.XGRID.size))                      # uninformative observation -> score 0
    lp = HO.gaussian_logp([3.0], [0.5])
    assert abs(HO.log_score(lp, L)[0]) < 1e-3
    L2 = HO.observation_loglik(_df().iloc[[0]])
    good = HO.log_score(HO.gaussian_logp([3.0], [0.3]), L2)[0]
    bad = HO.log_score(HO.gaussian_logp([1.0], [0.3]), L2)[0]
    assert good > bad + 5


def test_gaussian_baseline_recovers_truth_and_uses_zero_counts():
    rng = np.random.default_rng(0)
    n = 80
    xv = rng.normal(2.0, 0.6, n)
    V = 0.016
    N = rng.poisson(10 ** xv * V)
    df = pd.DataFrame(dict(sample=[str(i) for i in range(n)], cells=np.where(N > 0, N / V, 0.0),
                           cells_counted=N.astype(float), dataset="s6_2017"))
    L = HO.observation_loglik(df, v_zero=V)
    b, s, ok = HO.fit_gaussian_model(np.ones((n, 1)), L)
    assert ok and abs(b[0] - 2.0) < 0.2 and abs(s - 0.6) < 0.15
    # dropping zeros (old behaviour) biases the mean upwards
    pos = N > 0
    b2, _, _ = HO.fit_gaussian_model(np.ones((pos.sum(), 1)), L[pos])
    assert (N == 0).sum() > 3 and b2[0] > b[0]


def test_physics_fold_uses_training_site_only():
    nodes = np.round(np.arange(0.0, 6.0001, 0.05), 3)
    n, H = 6, 3
    logz = np.zeros((n, H))
    logz[:3, 0] = 5.0                                     # training samples prefer hyper 0
    logz[3:, 2] = 50.0                                    # test samples would prefer hyper 2 (must be ignored)
    marg = np.zeros((n, H, nodes.size))
    marg[:, :, 60] = 1.0                                  # all mass at log B = 3
    tab = dict(logz=logz, marg=marg, nodes=nodes, hyp=[(0.01, "a", 0, 1), (0.02, "b", 0, 1), (0.03, "c", 0, 1)],
               fwd=np.zeros((n, H)), bb=np.zeros((n, H, 2)))
    df = pd.DataFrame(dict(sample=list("abcdef"), cells=[1e3, 1e4, 1e2, 1e3, 1e3, 1e3],
                           cells_counted=[16, 160, 1.6, 16, 16, 16], dataset=["s6_2017"] * 6))
    L = HO.observation_loglik(df)
    train = np.array([1, 1, 1, 0, 0, 0], bool)
    r = HO.physics_fold(tab, train, ~train, L)
    assert r["h"] == 0
    assert r["tau"] > 0.5                                 # training obs spread 2 dex around a point prediction
    r2 = HO.physics_fold(dict(tab, logz=np.where(np.arange(n)[:, None] < 3, logz, 0.0)), train, ~train, L)
    assert r2["h"] == r["h"] and r2["tau"] == r["tau"]   # changing test-site evidence changes nothing


def test_fast_tables_equal_reference_grid_posterior(monkeypatch):
    import dataclasses
    import emulator as E
    import inversion as INV
    from priors import PriorConfig, prior_logpdfs
    rng = np.random.default_rng(1)
    axes = {"log_b": np.round(np.arange(1, 5.001, 0.5), 3), "f_n": np.linspace(0, 1, 3),
            "r_um": np.geomspace(500, 8000, 5), "dust_ppb": np.geomspace(3e4, 1e6, 3)}
    shp = tuple(len(v) for v in axes.values())
    g = np.meshgrid(*axes.values(), indexing="ij")
    base = 0.8 - 0.06 * g[0] - 0.02 * g[1] - 0.03 * np.log(g[2] / 500) - 0.01 * np.log(g[3] / 3e4)
    data = {"bands": np.stack([base * c for c in (1.0, 0.97, 0.94, 0.8)], -1), "bba": base * 0.9}
    em = E.Emulator(axes, data, {})
    monkeypatch.setattr(E, "build_emulator", lambda *a_, **k: em)
    monkeypatch.setattr(E.Emulator, "refine_log_b", lambda self, step: self)
    monkeypatch.setattr(HO.FV, "SIGMA_GRID", np.array([0.02, 0.04]))
    monkeypatch.setattr(HO.FV, "EB_MU_LNR", np.log([800.0, 3000.0]))
    monkeypatch.setattr(HO.FV, "EB_SD_LNR", np.array([0.5, 1.0]))
    df = pd.DataFrame(dict(sample=list("abcd"), cells=[1e3, 1e2, 0.0, 3e3], cells_counted=[16, 1.6, 0, 48],
                           dataset="s6_2017", sza=[45.2, 45.4, 47.0, 47.1], B2=0.6 + 0.02 * rng.normal(size=4)))
    for b, c in zip(("B3", "B4", "B8"), (0.97, 0.94, 0.8)):
        df[b] = df.B2 * c
    tab = HO.physics_tables(df, "tddft_D", None, None, 1, ".", verbose=False)
    pc0 = PriorConfig.for_density(690.0)
    R = df[HO.BANDS].to_numpy()
    for j in (0, 3, len(tab["hyp"]) - 1):
        s_, _, m_, d_ = tab["hyp"][j]
        lp, sk, mk, _ = prior_logpdfs(axes, 4, dataclasses.replace(pc0, mu_lnr=m_, sd_lnr=d_))
        ref = INV.GridPosterior(em, np.full(4, s_), derived=()).run(R, lp, sk, mk=mk, marginals=("log_b",))
        assert np.allclose(tab["logz"][:, j], ref["log_evidence"], atol=1e-8)
        assert np.allclose(tab["marg"][:, j], ref["marg_log_b"], atol=1e-10)
        # forward density = evidence with a delta prior at the node of the observed count
        y = np.log10(df.cells[0])
        lpd = dict(lp, log_b=np.where(np.isclose(axes["log_b"], np.round(y * 2) / 2), 0.0, -np.inf)[None, :])
        rf = INV.GridPosterior(em, np.full(4, s_), derived=()).run(R[:1], lpd, sk, mk=mk)
        assert np.isclose(tab["fwd"][0, j], rf["log_evidence"][0])
    assert np.isnan(tab["fwd"][2]).all()                       # zero count: no forward test


def test_requested_observation_sd_reaches_all_metrics():
    df = pd.DataFrame(dict(sample=["210805-S1", "210806-S2"], cells=[1000.0, 3000.0], cells_counted=[np.nan, np.nan],
                           dataset="sgris_2021"))
    lp = HO.gaussian_logp([3.0, 3.0], [0.05, 0.05])
    out = {}
    for sd in (0.01, 0.5):
        L = HO.observation_loglik(df, obs_sd_sgris=sd)
        out[sd] = HO.summary_stats(lp, df, L, obs_sd_sgris=sd)
    # sample 2 (log 3.48) is outside a narrow predictive but inside once a 0.5 dex observation error is used
    assert out[0.01]["covered"][1] == 0 and out[0.5]["covered"][1] == 1
    assert out[0.5]["crps"][1] < out[0.01]["crps"][1]
    assert not np.allclose(out[0.01]["log_score"], out[0.5]["log_score"])


def test_ridge_standardisation_is_fitted_inside_inner_folds(monkeypatch):
    rng = np.random.default_rng(0)
    n = 12
    df = pd.DataFrame(dict(sample=[f"{i}_7_SB{i}" for i in range(n)], cells=10 ** rng.uniform(2, 4, n),
                           cells_counted=50.0, dataset=["s6_2017"] * 6 + ["sgris_2021"] * 6))
    for b in HO.BANDS:
        df[b] = rng.uniform(0.3, 0.7, n)
    L = HO.observation_loglik(df)
    seen = []
    real = HO.fit_gaussian_model

    def spy(X, Lx, l2=0.0):
        if X.shape[1] > 3:
            seen.append(X[:, 1:].mean(axis=0))
        return real(X, Lx, l2)
    monkeypatch.setattr(HO, "fit_gaussian_model", spy)
    train = (df.dataset == "s6_2017").to_numpy()
    HO.baselines(df, train, ~train, L)
    # every ridge fit sees features standardised on exactly its own fitting rows (column means ~ 0)
    assert seen and all(np.allclose(m, 0.0, atol=1e-9) for m in seen)


def test_block_bootstrap_uses_days_and_flags_few_blocks():
    a = np.array([1.0, 1.0, 1.0, -1.0, -1.0, -1.0])
    b = np.zeros(6)
    r = HO.block_bootstrap(a, b, np.array(["d1"] * 3 + ["d2"] * 3), n=2000)
    assert r["n_blocks"] == 2 and r["few_blocks"] and r["lo"] <= -0.99 and r["hi"] >= 0.99
    df = pd.DataFrame(dict(sample=["13_7_SB5", "21_7_SB1", "210805-S1"], dataset=["s6_2017", "s6_2017", "sgris_2021"]))
    assert list(HO.sampling_block(df)) == ["13_7", "21_7", "210805"]
