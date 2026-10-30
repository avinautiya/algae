"""Regression for the small-thermal-exposure cancellation found by the independent audit (ice_network_verification).
The stress setting (a ~ 1e-13) is a numerical test, not a plausible heat-transfer closure; default runs are not claimed wrong.
Reference: 60-digit Decimal evaluation of the exact segment solution; legacy: the equilibrium-temperature form used before."""
import math
from decimal import Decimal, getcontext

import numpy as np
import pytest

from ice_network.constants import CP_W, RHO_W
from ice_network.thermal.channel_heat import heat_transfer_coeff, reach_heat, segment_heat

getcontext().prec = 60
D = Decimal


def exact(c, ua, theta, phi):
    c, ua, theta, phi = D(c), D(ua), D(theta), D(phi)
    if ua == 0:
        return theta + phi / c, D(0)
    a = ua / c
    e = (-a).exp()
    f = 1 - e
    return theta * e + phi / c * (f / a), c * theta * f + phi * (1 - f / a)


def legacy(c, ua, theta, phi):
    """The equilibrium-temperature form used by the first release (kept only to demonstrate the defect)."""
    hp_over_c = ua / c
    t_eq = phi / ua if ua > 0 else math.inf
    decay = math.exp(-hp_over_c)
    t_next = t_eq + (theta - t_eq) * decay
    return t_next, c * (theta - t_next) + phi


CASES = [(c, a, th, phi) for c in (2e3, 2e5) for a in (0.0, 1e-16, 1e-13, 1e-10, 1e-8, 1e-5, 1e-4, 1e-3, 0.1, 0.49, 0.51, 1.0, 5.0, 50.0)
         for th in (0.0, 0.3) for phi in (0.0, 1.0, 40.0)]


@pytest.mark.parametrize("c,a,theta,phi", CASES)
def test_segment_heat_matches_60_digit_reference_across_thermal_exposure(c, a, theta, phi):
    to, q = segment_heat(c, a * c, theta, phi)
    to_ref, q_ref = exact(c, a * c, theta, phi)
    scale = max(abs(float(D(c) * D(theta))), phi, 1e-30)           # watts-scale of the problem
    assert abs(q - float(q_ref)) <= 1e-13 * scale + 1e-15 * abs(float(q_ref)), (q, float(q_ref))
    assert abs(to - float(to_ref)) <= 1e-14 * max(abs(theta), phi / c, 1e-30)
    assert q >= 0.0 and to >= 0.0                                     # no clipping was applied; the exact values are nonnegative


def test_legacy_formulation_fails_at_the_stress_exposure_and_the_new_one_does_not():
    c, ua, theta, phi = 2e5, 1e-13 * 2e5, 0.3, 1.0
    q_ref = float(exact(c, ua, theta, phi)[1])
    assert q_ref > 0
    assert legacy(c, ua, theta, phi)[1] < 0 or abs(legacy(c, ua, theta, phi)[1] - q_ref) > 0.1 * q_ref   # defect reproduced
    assert segment_heat(c, ua, theta, phi)[1] == pytest.approx(q_ref, rel=1e-12)


def test_ordinary_conditions_are_unchanged_by_the_new_algebra():
    c, ua, theta, phi = 2e5, 0.3 * 2e5, 0.3, 40.0
    a = segment_heat(c, ua, theta, phi); b = legacy(c, ua, theta, phi)
    assert a[0] == pytest.approx(b[0], rel=1e-12) and a[1] == pytest.approx(b[1], rel=1e-9)


def test_nonuniform_geometry_unequal_segments_and_entrance_dissipation_against_chained_reference():
    q = 0.02
    h = np.array([0.30, 0.22, 0.15, 0.11, 0.09]); w = np.array([0.4, 0.45, 0.6, 0.5, 0.35])
    for dx in (3.0, 11.0):
        for phi_in, phi_seg in ((0.0, np.zeros(4)), (7.0, np.array([2.0, 5.0, 1.0, 9.0]))):
            for nu in (1.0, 1e-9, 1e-12):
                t_arr = np.zeros(5); qi = np.zeros(4)
                reach_heat(q, 0.25, 0.0, h, w, dx, phi_seg, phi_in, nu, t_arr, qi)
                theta = D(0.25); c = D(RHO_W * CP_W * q); q_ref = []
                for k in range(4):
                    wm, hm = 0.5 * (w[k] + w[k + 1]), 0.5 * (h[k] + h[k + 1])
                    ua = heat_transfer_coeff(q, wm, hm, nu) * (wm + 2 * hm) * dx
                    ph = phi_seg[k] + (phi_in if k == 0 else 0.0)          # entrance dissipation enters segment 0 only
                    theta, qk = exact(c, ua, theta, ph)
                    q_ref.append(float(qk))
                assert np.allclose(qi, q_ref, rtol=1e-9, atol=1e-13)
                assert t_arr[-1] == pytest.approx(float(theta), rel=1e-11, abs=1e-14)
                assert (qi >= 0).all()


def test_zero_conductance_limit_returns_all_dissipation_as_sensible_heat():
    to, q = segment_heat(2e5, 0.0, 0.3, 5.0)
    assert to == pytest.approx(0.3 + 5.0 / 2e5, rel=1e-15) and q == 0.0
