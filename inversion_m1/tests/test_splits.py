"""Split-leakage checks: synthetic cases and the real committed sample-group assignment."""

import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import splits as S  # noqa: E402

GROUPS_CSV = os.path.join(HERE, "..", "..", "records", "inversion_m1", "sample_groups.csv")
MANIFEST = os.path.join(HERE, "..", "..", "records", "inversion_m1", "data_manifest.json")


# ----------------------------------------------------------------------------- synthetic behaviour
def test_leave_one_group_out_covers_every_group_once_without_overlap():
    g = np.array(["a", "a", "b", "c", "c", "c"])
    seen = []
    for u, tr, te in S.leave_one_group_out(g):
        seen.append(u)
        assert set(g[te]) == {u} and u not in set(g[tr])
        assert len(set(tr) & set(te)) == 0 and len(tr) + len(te) == len(g)
    assert seen == ["a", "b", "c"]


def test_assert_disjoint_flags_leakage():
    S.assert_disjoint(["x", "y"], ["z"])
    with pytest.raises(AssertionError, match="leakage"):
        S.assert_disjoint(["x", "y"], ["y", "z"], "dev/test")


def test_train_only_blocks_held_out_reads():
    data = np.arange(10.0)
    v = S.TrainOnly(data, [0, 1, 2, 3])
    assert v[2] == 2.0 and v[np.array([0, 3])].tolist() == [0.0, 3.0]
    with pytest.raises(PermissionError):
        v[7]
    with pytest.raises(PermissionError):
        v[np.array([1, 8])]
    assert v.reads == {0, 2, 3}


def test_validate_roles_rejects_mixed_roles_duplicates_and_unknown():
    base = pd.DataFrame(dict(dataset=["d", "d", "d"], sample_id=["1", "2", "3"], site=["s"] * 3,
                             sampling_group=["g1", "g1", "g2"], role=["development", "development", "final_test_withheld"]))
    S.validate_roles(base)
    with pytest.raises(ValueError, match="mixed roles"):
        S.validate_roles(base.assign(role=["development", "final_test_withheld", "development"]))
    with pytest.raises(ValueError, match="duplicate"):
        S.validate_roles(pd.concat([base, base.iloc[[0]]]))
    with pytest.raises(ValueError, match="unknown roles"):
        S.validate_roles(base.assign(role=["development", "development", "oops"]))
    shared = pd.concat([base.iloc[[2]].assign(dataset="e", sample_id="9", role="development")])
    with pytest.raises(ValueError, match="different roles"):
        S.validate_roles(pd.concat([base, shared]))


def test_roles_hash_is_order_independent_and_sensitive_to_role():
    d = pd.DataFrame(dict(dataset=["d", "d"], sample_id=["1", "2"], sampling_group=["g", "g"], role=["development"] * 2))
    assert S.roles_hash(d) == S.roles_hash(d.iloc[::-1])
    assert S.roles_hash(d) != S.roles_hash(d.assign(role=["development", "final_test_withheld"]))


# --------------------------------------------------------------------------------------- real data
@pytest.fixture(scope="module")
def groups():
    return S.load_groups(GROUPS_CSV)


def test_real_assignment_is_structurally_valid(groups):
    S.validate_roles(groups)
    assert set(groups.role) <= S.ROLES_ALLOWED


def test_real_development_and_transfer_sets_are_disjoint_by_sample_and_group(groups):
    dev = groups[groups.role.isin(["development", "development_unlabelled"])]
    xfer = groups[groups.role == "transfer_check_reused"]
    ext = groups[groups.role == "external_method_reproduction"]
    S.assert_disjoint(dev.sampling_group, xfer.sampling_group, "development/transfer sampling groups")
    S.assert_disjoint(dev.sampling_group, ext.sampling_group, "development/external sampling groups")
    S.assert_disjoint(dev.site, xfer.site, "development/transfer sites")
    assert xfer.sampling_group.nunique() == 2 and dev[dev.role == "development"].sampling_group.nunique() == 9


def test_no_withheld_final_test_sample_is_present_in_the_working_tables(groups):
    # the withheld PROMBIO 2024 station-days must not appear in any table used for fitting or model selection
    assert (groups.role == "final_test_withheld").sum() == 0
    assert not groups.dataset.str.contains("PROMBIO").any()


def test_count_bearing_development_plots_are_exactly_the_colocated_ones(groups):
    dev = groups[groups.role == "development"]
    assert dev.has_counts.all() and dev.has_albedo.all() and dev.has_hcrf.all() and len(dev) == 47
    unl = groups[groups.role == "development_unlabelled"]
    assert (~unl.has_counts).all() and len(unl) == 39


def test_table_s4_rows_that_are_s6_plots_are_flagged_and_excluded_from_the_2021_set(groups):
    import json
    m = json.load(open(MANIFEST))
    sg = next(d for d in m["datasets"] if d["id"] == "SGreenland_2021_hcrf_counts")
    s6_rows = set(sg["table_s4_rows_that_are_S6_2017_plots"])
    assert len(s6_rows) == 6
    x = groups[groups.dataset == "SGreenland_2021_hcrf_counts"]
    S.assert_disjoint(s6_rows, set(x.sample_id), "Table S4 S6-2017 rows / 2021 evaluation spectra")
    # ... but those 6 rows ARE S6 plots, i.e. they share identity with the development set (naming differs: '-' vs '_')
    dev_ids = {i.replace("_", "-") for i in groups[groups.role == "development"].sample_id}
    assert s6_rows <= dev_ids


def test_hyperparameter_selection_cannot_read_the_held_out_day(groups):
    dev = groups[groups.role == "development"].reset_index(drop=True)
    for held, tr, te in S.leave_one_group_out(dev.sampling_group):
        view = S.TrainOnly(np.arange(len(dev)), tr)
        _ = [view[int(i)] for i in tr]                         # training reads are allowed
        assert view.reads == set(int(i) for i in tr)
        with pytest.raises(PermissionError):
            view[int(te[0])]
