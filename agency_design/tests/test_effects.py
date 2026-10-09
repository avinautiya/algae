import unittest
import numpy as np
from agency_design.effects import effects_from_albedo, aggregate_pixels, rank_effects


class EffectTests(unittest.TestCase):
    def gates(self, regions):
        return {r: dict(observation_qc_pass=True, in_domain=True,
                        retrieval_identifiable=True, forward_validation_pass=True) for r in regions}

    def rank(self, q, quality=None):
        return rank_effects(q, ["a", "b"], list(range(len(q))), quantity="absorbed_sw", units="W m-2",
                            minimum_effect=1, practical_difference=2,
                            quality=self.gates(["a", "b"]) if quality is None else quality)

    def test_spectral_effect_signed(self):
        r = effects_from_albedo(np.array([[[.6, .3]]]), np.array([[[.4, .4]]]), np.array([[[100, 50]]]))
        np.testing.assert_allclose(r, [[15]])

    def test_brightening_not_clipped(self):
        r = effects_from_albedo(np.array([[[.3]]]), np.array([[[.4]]]), np.array([[[100]]]))
        self.assertLess(r[0, 0], 0)

    def test_joint_aggregation_preserves_shared_error(self):
        r = aggregate_pixels([[10, 10], [20, 20]], [2, 3], quantity="absorbed_sw_w_m2")
        np.testing.assert_equal(r["total"], [50, 100])
        np.testing.assert_equal(r["mean_per_area"], [10, 20])

    def test_resolved_difference_with_shared_uncertainty(self):
        r = self.rank([[10, 5], [100, 95], [1000, 995]])
        self.assertEqual(r["pairwise"][0]["higher_effect_region"], "a")
        np.testing.assert_allclose(r["pairwise"][0]["difference_interval95"], [5, 5])

    def test_unresolved_not_forced_total_order(self):
        r = self.rank([[10, 9], [9, 10]])
        self.assertEqual(r["pairwise"][0]["ranking"], "unresolved")
        self.assertIsNone(r["total_order"])

    def test_failed_identifiability_overrides_large_effect(self):
        gates = self.gates(["a", "b"])
        gates["a"]["retrieval_identifiable"] = False
        r = self.rank([[100, 10], [110, 11]], gates)
        self.assertEqual(r["regions"][0]["status"], "unsupported")
        self.assertIsNone(r["pairwise"][0]["higher_effect_region"])

    def test_missing_gate_rejected(self):
        gates = self.gates(["a", "b"])
        del gates["a"]["forward_validation_pass"]
        with self.assertRaises(ValueError):
            self.rank([[1, 2], [3, 4]], gates)

    def test_missing_draws_rejected(self):
        with self.assertRaises(ValueError):
            self.rank([[np.nan, 3], [2, 4]])


if __name__ == "__main__":
    unittest.main()
