"""H2/H4 scoring: station-year cluster contrasts, the < 5 block rule, verdict logic (synthetic rows)."""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))

import h2h4_station as H  # noqa: E402

MODELS = ["M0", "M1", "M1b", "M2", "M3", "M3alt_C", "M3alt_iid"]


def _df(years, m3_better_by=0.03, n_per=6, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for st in ("KAN_L", "KAN_M"):
        for y in years:
            for _ in range(n_per):
                obs = rng.uniform(0.3, 0.6)
                sw = 700.0
                r = dict(station=st, year=y, albedo_obs=obs, sw_down=sw, absorbed_obs=sw * (1 - obs),
                         excluded_node=False, node_offset=0.0)
                for m in MODELS:
                    e = rng.normal(0, 0.005) + (0.0 if m == "M3" else m3_better_by)
                    r[f"alpha_{m}"] = obs + e
                rows.append(r)
    df = pd.DataFrame(rows)
    for m in MODELS:
        df[f"abserr_absorbed_{m}"] = (df.sw_down * (1 - df[f"alpha_{m}"]) - df.absorbed_obs).abs()
        df[f"abserr_albedo_{m}"] = (df[f"alpha_{m}"] - df.albedo_obs).abs()
        df[f"err_albedo_{m}"] = df[f"alpha_{m}"] - df.albedo_obs
    return df


def test_supported_with_enough_blocks():
    S, C, v = H.evaluate(_df((2017, 2018, 2020)))                  # 6 primary blocks
    assert v["H4"].startswith("SUPPORTED") and v["H2_vs_M2"].startswith("SUPPORTED")
    assert "H1" in v["H2_vs_M2"]                                     # H2 needs H1 not contradicted


def test_too_few_blocks_gives_no_interval():
    S, C, v = H.evaluate(_df((2017, 2018)))                         # 4 blocks
    assert v["H4"].startswith("INSUFFICIENT EVIDENCE")
    row = C[(C.set == "primary") & (C.hypothesis == "H4")].iloc[0]
    assert np.isnan(row.lo) and row.n_blocks == 4


def test_below_minimum_is_not_supported():
    S, C, v = H.evaluate(_df((2017, 2018, 2020), m3_better_by=0.004))   # real but < 0.01 albedo
    assert v["H4"].startswith("NOT SUPPORTED")


def test_2019_kept_separate():
    S, C, v = H.evaluate(_df((2017, 2018, 2019, 2020)))
    assert set(C.set) >= {"primary", "2019_separate"}
    assert S[S.set == "primary"].n.iloc[0] == 2 * 3 * 6
