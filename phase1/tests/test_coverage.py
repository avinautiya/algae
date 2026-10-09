"""Regression: root coverage is judged over the window consumed downstream (260-750 nm, 265-600 nm
normalisation), not at the 300 nm plotting edge."""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

import spectra  # noqa: E402


def test_window_matches_calibration():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "phase2"))
    import empirical_data as ED
    import tddft_calibration as tc
    assert spectra.DOWNSTREAM_WINDOW_NM[0] <= tc.MAC_RANGE_NM[0] and spectra.DOWNSTREAM_WINDOW_NM[1] >= tc.MAC_RANGE_NM[1]
    assert tuple(spectra.CAL_NORM_WINDOW_NM) == tuple(ED.FE_NORM_RANGE_NM)


def test_root_inside_window_is_flagged_root_far_above_is_not():
    e_in = np.array([2.5, 3.5, 4.6])            # top root at 270 nm: inside the 265-600 nm window
    e_far = np.array([2.5, 3.5, 9.0])           # top root at 138 nm
    f = np.array([0.1, 0.2, 0.5])
    a = spectra.root_count_sensitivity(e_in, f, 300.0, 0.3, drop=1)
    b = spectra.root_count_sensitivity(e_far, f, 300.0, 0.3, drop=1)
    assert a["rel_change_norm_integral"] > 0.05 and a["rel_change_at_short_edge"] > 0.5
    assert b["max_rel_change_in_window"] < 1e-6
    # the old criterion (highest root + 2 FWHM reaches 300 nm) passes the first case, the new one does not
    assert e_in.max() + 2 * 0.3 >= spectra.HC_EV_NM / 300.0
