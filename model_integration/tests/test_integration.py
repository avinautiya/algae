"""Analytic tests, NOT empirical evidence of glacier-model improvement."""
import unittest
import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from model_integration.coupling import (
    spectral_budget, couple_host_albedo, mix_subpixels, integrate_intervals,
    melt_diagnostic, paired_melt_diagnostic, LATENT_HEAT_J_KG,
)
from model_integration.benchmark import score, audit_split, REQUIRED_MODELS
from model_integration.biosnicar_adapter import AlgaeOpticsAdapter
from model_integration.export import export_spectra, load_spectra
from model_integration.benchmark import main as benchmark_main


class CouplingTests(unittest.TestCase):
    def test_spectral_energy(self):
        out = spectral_budget([.6, .4], [.4, .3], [200, 100])
        self.assertAlmostEqual(out["absorbed_without_w_m2"], 140)
        self.assertAlmostEqual(out["absorbed_with_w_m2"], 190)
        self.assertAlmostEqual(out["delta_absorbed_w_m2"], 50)

    def test_zero_and_brightening(self):
        self.assertEqual(spectral_budget([.5], [.5], [0])["delta_absorbed_w_m2"], 0)
        self.assertLess(spectral_budget([.3], [.4], [100])["delta_absorbed_w_m2"], 0)

    def test_invalid_spectrum(self):
        for a, f in (([np.nan], [1]), ([1.1], [1]), ([.5], [-1]), ([.5], [1, 2])):
            with self.assertRaises(ValueError):
                spectral_budget(a, a, f)

    def test_double_count_guard(self):
        for treatment in ("explicit", "implicit_observed", "unknown"):
            with self.assertRaises(ValueError):
                couple_host_albedo([.5], [.6], [.4], mode="anomaly", host_algae_treatment=treatment)
        np.testing.assert_allclose(couple_host_albedo([.5], [.6], [.4], mode="anomaly",
                                                     host_algae_treatment="absent"), [.3])
        np.testing.assert_allclose(couple_host_albedo([.5], [.6], [.4], mode="replace",
                                                     host_algae_treatment="implicit_observed"), [.4])

    def test_no_silent_clipping(self):
        with self.assertRaises(ValueError):
            couple_host_albedo([.1], [.9], [.1], mode="anomaly", host_algae_treatment="absent")

    def test_subpixel_mixture(self):
        np.testing.assert_allclose(mix_subpixels([[.8, .6], [.2, .4]], [.25, .75]), [.35, .45])
        with self.assertRaises(ValueError):
            mix_subpixels([[.8], [.2]], [.2, .2])

    def test_time_energy(self):
        self.assertEqual(integrate_intervals([100, 200], [3600, 1800]), 720000)
        with self.assertRaises(ValueError):
            integrate_intervals([100], [0])

    def test_cold_content_and_energy_closure(self):
        out = melt_diagnostic([100, 0, 100], [0, -50, 0], [10, 10, 10],
                              initial_cold_content_j_m2=600,
                              surface_at_melting_point=np.array([True, True, True]))
        np.testing.assert_allclose(out["potential_melt_mm_we"] * LATENT_HEAT_J_KG, [400, 0, 500])
        np.testing.assert_allclose(out["cold_content_j_m2"], [0, 500, 0])
        self.assertAlmostEqual(1500 + out["cold_content_j_m2"][-1] - 600,
                               out["potential_melt_mm_we"].sum() * LATENT_HEAT_J_KG)

    def test_not_at_melting_point(self):
        out = melt_diagnostic([100], [0], [10], initial_cold_content_j_m2=0,
                              surface_at_melting_point=np.array([False]))
        self.assertEqual(out["potential_melt_mm_we"][0], 0)
        self.assertEqual(out["unallocated_energy_j_m2"][0], 1000)

    def test_paired_melt(self):
        out = paired_melt_diagnostic([100], [150], [-120], [3600],
                                    initial_cold_content_j_m2=0,
                                    surface_at_melting_point=np.array([True]))
        self.assertAlmostEqual(out["delta_potential_melt_mm_we"][0], 30 * 3600 / LATENT_HEAT_J_KG)


class BenchmarkTests(unittest.TestCase):
    def test_group_leakage(self):
        with self.assertRaises(ValueError):
            audit_split(["train"], ["test"], ["S6-2017"], ["S6-2017"])
        with self.assertRaises(ValueError):
            audit_split([], ["test", "test"], [], ["a"])
        audit_split(["a"], ["b"], ["A"], ["B"])

    def test_paired_bootstrap_reproducible(self):
        pred = {k: [1, 1, 1, 1] for k in REQUIRED_MODELS}
        pred["pigment_cell"] = [0, 0, 0, 0]
        args = ([0, 0, 0, 0], pred, ["a", "a", "b", "b"])
        a = score(*args, bootstrap=100)
        self.assertEqual(a, score(*args, bootstrap=100))
        self.assertEqual(a["n_groups"], 2)
        self.assertEqual(a["rmse_improvement_vs_baseline"]["established_algae"]["ci95"], [1, 1])

    def test_single_group_no_interval(self):
        out = score([0], {k: [0] for k in REQUIRED_MODELS}, ["a"], bootstrap=100)
        self.assertIsNone(out["rmse_improvement_vs_baseline"]["no_algae"]["ci95"])

    def test_missing_predictions_fail(self):
        pred = {k: [0] for k in REQUIRED_MODELS}
        pred["pigment_cell"] = [np.nan]
        with self.assertRaises(ValueError):
            score([0], pred, ["a"], bootstrap=100)


class AdapterTests(unittest.TestCase):
    def factory(self, cfg, *_args):
        class Runner:
            direct = 1
            _ill_cache = {}
            band = np.array([True, True])
            wvl_um = np.array([.5, 1.0])
            def run_multi(self, spec, sza, imps):
                B = sum(c for name, c in imps if name != "dust")
                return np.array([.6 - .01 * B, .4]), np.array([2., 1.]), 0
        return SimpleNamespace(cfg=cfg, runner=Runner(), dust="dust",
                               imps={"nordenskioeldii": "n", "alaskanum": "a"},
                               spec=lambda r: SimpleNamespace(split=False))

    def test_adapter_matches_budget_and_zero_algae(self):
        adapter = AlgaeOpticsAdapter(SimpleNamespace(model="ours"), builder_factory=self.factory)
        kw = dict(f_n=.5, r_um=1000, dust_ppb=10, sza_deg=45, sw_down_w_m2=300)
        out = adapter.evaluate(cells_ml_meltwater=10, **kw)
        np.testing.assert_allclose(out["flux_bin_w_m2"], [200, 100])
        self.assertAlmostEqual(out["delta_absorbed_w_m2"], 20)
        self.assertAlmostEqual(adapter.evaluate(cells_ml_meltwater=0, **kw)["delta_absorbed_w_m2"], 0)

    def test_night_rejected(self):
        adapter = AlgaeOpticsAdapter(SimpleNamespace(model="ours"), builder_factory=self.factory)
        with self.assertRaises(ValueError):
            adapter.evaluate(cells_ml_meltwater=0, f_n=.5, r_um=1000,
                             dust_ppb=0, sza_deg=95, sw_down_w_m2=10)


class ExportAndCliTests(unittest.TestCase):
    def test_roundtrip_checksum_and_overwrite(self):
        provenance = dict(source_commit="test-only", input_sha256={}, biosnicar_commit="mock",
                          state_coordinates={}, leading_dimensions=[], illumination="mock",
                          calibration_draw_ids=[], validation_status="analytic test only")
        with tempfile.TemporaryDirectory() as directory:
            dst = Path(directory) / "spectra.npz"
            export_spectra(dst, [500, 1000], [.6, .4], [.4, .3], [200, 100], provenance=provenance)
            arrays, metadata = load_spectra(dst)
            np.testing.assert_equal(arrays["wavelength_nm"], [500, 1000])
            self.assertEqual(metadata["shape"], [2])
            with self.assertRaises(FileExistsError):
                export_spectra(dst, [500, 1000], [.6, .4], [.4, .3], [200, 100], provenance=provenance)
            with dst.open("ab") as handle:
                handle.write(b"tampered")
            with self.assertRaises(ValueError):
                load_spectra(dst)

    def test_cli_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            csv = root / "predictions.csv"
            csv.write_text("sample_id,group_id,observed," + ",".join(REQUIRED_MODELS) +
                           "\nx,A,0,1,1,1,0\ny,B,0,1,1,1,0\n")
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps(dict(quantity="test-only", units="1",
                observation_source="synthetic analytic test, not field evidence", training_ids=["z"],
                training_groups=["C"], model_provenance={k: "mock" for k in REQUIRED_MODELS},
                calibration_frozen_before_test=True)))
            output = root / "report.json"
            benchmark_main([str(csv), str(manifest), "--output", str(output)])
            report = json.loads(output.read_text())
            self.assertEqual(report["n_groups"], 2)
            self.assertEqual(report["metrics"]["pigment_cell"]["rmse"], 0)


if __name__ == "__main__":
    unittest.main()
