import unittest
import numpy as np

from comparison_study.study import build, tune, paired_rmse, metrics, sampling_day


class StudyTests(unittest.TestCase):
    def fixture(self):
        rng = np.random.default_rng(3)
        X = rng.normal(size=(12, 2))
        return X, (X[:, :1] * .5), np.repeat(["a", "b", "c"], 4)

    def test_fold_local_scaling(self):
        X, y, g = self.fixture()
        fit = build("ridge", {"alpha": 1}).fit(X[g != "a"], y[g != "a"])
        np.testing.assert_allclose(fit.regressor_.steps[0][1].mean_, X[g != "a"].mean(0))
        np.testing.assert_allclose(fit.transformer_.mean_, y[g != "a"].mean(0))

    def test_tuning_has_no_test_input(self):
        X, y, g = self.fixture()
        fit, diag = tune(X, y, g, "ridge")
        self.assertEqual(diag["training_days"], ["a", "b", "c"])
        self.assertEqual(len(diag["candidates"]), 4)
        self.assertTrue(np.isfinite(fit.predict(X)).all())

    def test_all_models_multioutput(self):
        X, y, g = self.fixture()
        yy = np.column_stack([y, 2 * y])
        for model, params in (("mean", {}), ("ridge", {"alpha": 1}),
                              ("random_forest", {"min_samples_leaf": 2}),
                              ("gaussian_process", {"length_scale": 2, "noise": .2})):
            fit = build(model, params).fit(X, yy)
            self.assertEqual(fit.predict(X).shape, yy.shape)
            self.assertTrue(np.isfinite(fit.predict(X)).all())

    def test_single_group_rejected(self):
        X, y, _ = self.fixture()
        with self.assertRaises(ValueError):
            tune(X, y, np.repeat("a", len(X)), "ridge")

    def test_nonfinite_rejected(self):
        X, y, g = self.fixture()
        X[0, 0] = np.nan
        with self.assertRaises(ValueError):
            tune(X, y, g, "ridge")

    def test_few_blocks_no_interval(self):
        r = paired_rmse(np.zeros((4, 1)), np.zeros((4, 1)), np.ones((4, 1)), ["a", "a", "b", "b"])
        self.assertIsNone(r["ci95"])
        self.assertEqual(r["rmse_improvement"], 1)

    def test_block_interval_and_reproducibility(self):
        args = (np.zeros((5, 1)), np.zeros((5, 1)), np.ones((5, 1)), list("abcde"))
        a = paired_rmse(*args, n=100)
        self.assertEqual(a, paired_rmse(*args, n=100))
        self.assertEqual(a["ci95"], [1, 1])

    def test_exact_metrics(self):
        m = metrics(np.array([[0, 0]]), np.array([[1, -1]]))
        self.assertEqual(m["rmse"], 1)
        self.assertEqual(m["bias"], 0)
        self.assertEqual(m["per_output_bias"], [1, -1])

    def test_day_parsing(self):
        self.assertEqual(sampling_day("s6_2017", "20_7_SB4"), "s6_2017:2017-07-20")
        self.assertEqual(sampling_day("sgris_2021", "210805-S1"), "sgris_2021:2021-08-05")


if __name__ == "__main__":
    unittest.main()
