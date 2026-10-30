"""
Group-aware splits and leakage guards for the inversion study.

The unit of independence is the SAMPLING GROUP (a sampling day at one site), never the individual spectrum:
plots measured on the same day at the same site share illumination, ice state and operator effects. Hyper-
parameters (error scale, discrepancy) may only be chosen from TRAINING groups; a test group must never reach a
function that selects or tunes anything.
"""

from __future__ import annotations

import hashlib
import json
from typing import Iterable, Iterator

import numpy as np
import pandas as pd

ROLES_ALLOWED = {"development", "development_unlabelled", "transfer_check_reused", "external_method_reproduction",
                 "final_test_withheld"}


def load_groups(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    need = {"dataset", "sample_id", "site", "sampling_group", "role"}
    if not need <= set(df.columns):
        raise ValueError(f"sample-group table lacks columns {sorted(need - set(df.columns))}")
    return df


def validate_roles(df: pd.DataFrame) -> None:
    """Structural checks that make a split auditable."""
    bad = set(df.role) - ROLES_ALLOWED
    if bad:
        raise ValueError(f"unknown roles: {sorted(bad)}")
    dup = df[df.duplicated(["dataset", "sample_id"], keep=False)]
    if len(dup):
        raise ValueError(f"duplicate (dataset, sample_id): {dup[['dataset', 'sample_id']].head().values.tolist()}")
    # a sampling group must have exactly one role within a dataset, except development/unlabelled which share days
    g = df.assign(r=df.role.replace({"development_unlabelled": "development"})).groupby(["dataset", "sampling_group"]).r.nunique()
    if (g > 1).any():
        raise ValueError(f"sampling groups with mixed roles: {g[g > 1].index.tolist()}")
    # no sampling group may appear under two different datasets with different roles (a site shared across roles)
    s = df.assign(r=df.role.replace({"development_unlabelled": "development"})).groupby("sampling_group").r.nunique()
    if (s > 1).any():
        raise ValueError(f"a sampling group appears with different roles: {s[s > 1].index.tolist()}")


def assert_disjoint(a: Iterable, b: Iterable, what: str = "sets") -> None:
    inter = set(a) & set(b)
    if inter:
        raise AssertionError(f"leakage: {what} share {len(inter)} item(s), e.g. {sorted(inter)[:5]}")


def leave_one_group_out(groups) -> Iterator[tuple]:
    """Yield (held_out_group, train_index, test_index). Every group is held out exactly once; train and test
    never share a group."""
    g = np.asarray(groups)
    for u in sorted(set(g.tolist())):
        test = np.flatnonzero(g == u)
        train = np.flatnonzero(g != u)
        assert_disjoint(g[train], g[test], f"train/test groups (held out {u})")
        yield u, train, test


def roles_hash(df: pd.DataFrame) -> str:
    """Content hash of the frozen assignment (dataset, id, group, role), independent of row order."""
    rows = sorted(map(tuple, df[["dataset", "sample_id", "sampling_group", "role"]].astype(str).values.tolist()))
    return hashlib.sha256(json.dumps(rows).encode()).hexdigest()


class TrainOnly:
    """Wrap an array/table so that only the allowed row indices can be read. Hyper-parameter selection is given
    a TrainOnly view; reading a held-out row raises."""

    def __init__(self, data, allowed_index):
        self._d = data
        self._allowed = set(int(i) for i in allowed_index)
        self.reads: set[int] = set()

    def __getitem__(self, i):
        idx = np.atleast_1d(np.arange(len(self._d))[i]) if not isinstance(i, (int, np.integer)) else [int(i)]
        bad = [int(j) for j in idx if int(j) not in self._allowed]
        if bad:
            raise PermissionError(f"read of held-out row(s) {bad[:5]} during hyper-parameter selection")
        self.reads.update(int(j) for j in idx)
        return self._d[i]

    def __len__(self):
        return len(self._d)
