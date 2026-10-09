"""Host-coupling modes: declared algae treatment, double-counting refusals, Mode B flags, reference SEB."""
import os
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from model_integration import host_modes as HM  # noqa: E402


def _forcing(n=72):
    t = pd.date_range("2019-07-01", periods=n, freq="h")
    h = np.arange(n)
    return pd.DataFrame(dict(time=t, sw_down=np.clip(700 * np.sin(2 * np.pi * (h - 6) / 24), 0, None),
                             dlr=np.full(n, 290.0), T_a=274.0 + 2 * np.sin(2 * np.pi * (h - 9) / 24),
                             q_a=np.full(n, 0.004), wspd_u=np.full(n, 5.0), p_u=np.full(n, 870.0),
                             z_meas=np.full(n, 2.5)))


class TestContracts(unittest.TestCase):
    def test_treatment_must_be_declared(self):
        with self.assertRaises(ValueError):
            HM.HostContract("h", "unknown", "prognostic")
        with self.assertRaises(ValueError):           # observed albedo cannot be 'absent' of algae
            HM.HostContract("h", "absent", "observed")

    def test_mode_a_refuses_double_counting(self):
        obs = HM.HostContract("sat", "implicit_observed", "observed")
        with self.assertRaises(ValueError):
            HM.check_mode(obs, "A", "anomaly")
        expl = HM.HostContract("h", "explicit", "prognostic")
        for m in ("anomaly", "replace"):
            with self.assertRaises(ValueError):
                HM.check_mode(expl, "A", m)
        off = HM.HostContract("h", "explicit", "prognostic", own_algae_scheme_disabled=True)
        with self.assertRaises(ValueError):
            HM.check_mode(off, "A", "anomaly")
        HM.check_mode(off, "A", "replace")

    def test_mode_a_anomaly_on_algae_free_host(self):
        c = HM.HostContract("h", "absent", "prognostic")
        a, meta = HM.mode_a_albedo(c, np.array([0.5, 0.45]), np.array([0.55, 0.5]), np.array([0.5, 0.48]), "anomaly")
        np.testing.assert_allclose(a, [0.45, 0.43])
        self.assertEqual(meta["mode"], "A")

    def test_mode_b_keeps_observation_and_sign(self):
        c = HM.HostContract("s2", "implicit_observed", "observed")
        r = HM.mode_b_attribution(c, np.array([0.4, 0.4, 0.98]), np.array([0.05, -0.01, 0.05]))
        np.testing.assert_allclose(r["albedo"], [0.4, 0.4, 0.98])
        np.testing.assert_allclose(r["counterfactual_no_algae"][:2], [0.45, 0.39])   # brightening kept
        self.assertTrue(r["flag_counterfactual_out_of_range"][2])
        self.assertEqual(r["validation_status"], HM.NOT_INDEPENDENT)
        with self.assertRaises(ValueError):
            HM.check_mode(HM.HostContract("h", "absent", "prognostic"), "B")


class TestReferenceSEB(unittest.TestCase):
    def test_paired_runs_label_quantities(self):
        f = _forcing()
        c = HM.HostContract("reference SEB", "absent", "prognostic")
        r = HM.reference_seb_paired(c, f, np.full(len(f), 0.40), np.full(len(f), 0.45), mode="A")
        self.assertGreater(r["modelled_increment_mwe"], 0)
        self.assertGreater(r["potential_melt_mwe"], 0)
        self.assertEqual(r["runoff"], "not modelled")
        b = HM.reference_seb_paired(HM.HostContract("obs", "implicit_observed", "observed"), f,
                                    np.full(len(f), 0.40), np.full(len(f), 0.45), mode="B")
        self.assertAlmostEqual(b["modelled_increment_mwe"], r["modelled_increment_mwe"])
        self.assertEqual(b["validation_status"], HM.NOT_INDEPENDENT)


if __name__ == "__main__":
    unittest.main()
