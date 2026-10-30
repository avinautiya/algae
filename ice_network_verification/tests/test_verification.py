import copy
import csv
import json
import math
from pathlib import Path
import tempfile
import unittest

from ice_network_verification.__main__ import evaluate, generate, read_csv
from ice_network_verification.checks import (InvalidExport, audit_energy, audit_mass,
                                             compare_rows, observed_order)
from ice_network_verification.references import analytical_cases
from ice_network_verification.thermal import temperate_segment


def cases():
    return {case.name: case for case in analytical_cases()}


class AnalyticalReferences(unittest.TestCase):
    def test_linear_reservoir_half_life_and_integral(self):
        case = cases()["linear_reservoir"]
        row = case.solution(500 * math.log(2))
        self.assertAlmostEqual(row["storage_m3"], 5)
        self.assertAlmostEqual(row["discharge_m3_s"], 0.01)
        self.assertAlmostEqual(row["outflow_m3"], 5)

    def test_drying_has_known_time_and_no_negative_storage(self):
        case = cases()["finite_drainage"]
        self.assertEqual(case.parameters["emptying_time_s"], 200)
        self.assertAlmostEqual(case.solution(100)["level_m"], 0.25)
        for t in (200, 201, 1e6):
            self.assertEqual(case.solution(t)["storage_m3"], 0)
            self.assertEqual(case.solution(t)["discharge_m3_s"], 0)
            self.assertEqual(case.solution(t)["outflow_m3"], 4)

    def test_unequal_outlets_have_correct_ratio_not_hardcoded_symmetry(self):
        row = cases()["unequal_outlets"].solution(100)
        self.assertAlmostEqual(row["outlet_2_m3_s"] / row["outlet_1_m3_s"], 3)
        self.assertAlmostEqual(row["outlet_2_volume_m3"] / row["outlet_1_volume_m3"], 3)
        self.assertAlmostEqual(row["storage_m3"] + row["outlet_1_volume_m3"]
                               + row["outlet_2_volume_m3"], 6.4)

    def test_equal_outlets_are_exactly_symmetric(self):
        for row in cases()["symmetric_outlets"].rows():
            self.assertEqual(row["outlet_1_m3_s"], row["outlet_2_m3_s"])

    def test_sill_time_before_any_overtopping(self):
        case = cases()["pre_overtopping"]
        self.assertEqual(case.parameters["first_sill_contact_s"], 1000)
        self.assertAlmostEqual(case.solution(1000)["level_m"], 0.5)
        self.assertAlmostEqual(case.solution(1000)["inflow_volume_m3"], 2)

    def test_thermal_front_travel_time_and_attenuation(self):
        case = cases()["thermal_front"]
        self.assertEqual(case.solution(24.9)["temperature_c"], 0)
        self.assertAlmostEqual(case.solution(25)["temperature_c"], 2 * math.exp(-0.25))

    def test_heat_to_melt_distinguishes_solid_and_liquid_volume(self):
        row = cases()["heat_to_melt"].solution(3600)
        self.assertAlmostEqual(row["melt_mass_kg"] * 334000, row["heat_to_ice_j"])
        self.assertAlmostEqual(row["ice_volume_lost_m3"] * 917, row["melt_mass_kg"])
        self.assertAlmostEqual(row["meltwater_volume_m3"] * 1000, row["melt_mass_kg"])
        self.assertGreater(row["ice_volume_lost_m3"], row["meltwater_volume_m3"])
        self.assertAlmostEqual(row["outlet_temperature_c"], 2 * math.exp(-200 / 4180))


class ExportChecks(unittest.TestCase):
    def test_reference_identity_only_checks_harness(self):
        for case in analytical_cases():
            self.assertTrue(compare_rows(case.rows(), case.rows())["passed"])

    def test_wrong_ice_density_is_detected(self):
        reference = cases()["heat_to_melt"].rows()
        wrong = copy.deepcopy(reference)
        for row in wrong:
            row["ice_volume_lost_m3"] = row["meltwater_volume_m3"]
        result = compare_rows(reference, wrong)
        self.assertFalse(result["passed"])
        self.assertFalse(result["fields"]["ice_volume_lost_m3"]["passed"])

    def test_wrong_split_and_zero_discharge_are_detected(self):
        reference = cases()["unequal_outlets"].rows()
        wrong = copy.deepcopy(reference)
        for row in wrong:
            row["outlet_2_m3_s"] = row["outlet_1_m3_s"]
        self.assertFalse(compare_rows(reference, wrong)["passed"])
        wrong = copy.deepcopy(reference)
        for row in wrong:
            row["outlet_1_m3_s"] = 0
        self.assertFalse(compare_rows(reference, wrong)["passed"])

    def test_nonfinite_missing_negative_and_misaligned_fail_closed(self):
        reference = cases()["linear_reservoir"].rows()
        for value in (math.nan, math.inf, -0.1):
            wrong = copy.deepcopy(reference)
            wrong[1]["storage_m3"] = value
            with self.assertRaises(InvalidExport):
                compare_rows(reference, wrong)
        wrong = copy.deepcopy(reference)
        del wrong[1]["storage_m3"]
        with self.assertRaises(InvalidExport):
            compare_rows(reference, wrong)
        wrong = copy.deepcopy(reference)
        wrong[1]["time_s"] += 0.01
        with self.assertRaises(InvalidExport):
            compare_rows(reference, wrong)
        with self.assertRaises(InvalidExport):
            compare_rows(reference, reference[:-1])

    def test_duplicate_time_is_rejected(self):
        reference = cases()["linear_reservoir"].rows()
        wrong = copy.deepcopy(reference)
        wrong[1]["time_s"] = wrong[0]["time_s"]
        with self.assertRaises(InvalidExport):
            compare_rows(reference, wrong)

    def test_invalid_tolerances_rejected(self):
        ref = cases()["linear_reservoir"].rows()
        for rtol in (-1, math.inf, math.nan):
            with self.assertRaises(ValueError):
                compare_rows(ref, ref, rtol=rtol)


class BudgetsAndConvergence(unittest.TestCase):
    def test_nonzero_mass_flux_and_ice_conversion(self):
        rows = [dict(time_s=0, liquid_mass_kg=10, mass_in_kg=0,
                     mass_out_kg=0, melt_mass_kg=0),
                dict(time_s=1, liquid_mass_kg=12, mass_in_kg=3,
                     mass_out_kg=2, melt_mass_kg=1)]
        self.assertTrue(audit_mass(rows)["passed"])
        self.assertTrue(audit_mass(rows)["nonzero_inventory_change"])
        rows[1]["liquid_mass_kg"] += 0.083
        self.assertFalse(audit_mass(rows)["passed"])

    def test_no_flow_budget_is_labelled_degenerate(self):
        rows = [dict(time_s=t, liquid_mass_kg=10, mass_in_kg=0,
                     mass_out_kg=0, melt_mass_kg=0) for t in (0, 1)]
        self.assertTrue(audit_mass(rows)["passed"])
        self.assertFalse(audit_mass(rows)["nonzero_inventory_change"])

    def test_energy_audit_detects_double_counted_dissipation(self):
        rows = [dict(time_s=t, sensible_energy_j=0, latent_energy_j=15 * t,
                     boundary_heat_j=0, advected_energy_in_j=20 * t,
                     advected_energy_out_j=10 * t, dissipated_to_heat_j=5 * t)
                for t in (0, 1)]
        self.assertTrue(audit_energy(rows)["passed"])
        rows[1]["boundary_heat_j"] = 5
        self.assertFalse(audit_energy(rows)["passed"])

    def test_missing_ledger_field_rejected(self):
        with self.assertRaises(InvalidExport):
            audit_mass([{"time_s": 0}, {"time_s": 1}])

    def test_first_order_euler_against_independent_exponential(self):
        errors = []
        for n in (20, 40, 80, 160):
            value = 1.0
            dt = 1 / n
            for _ in range(n):
                value -= dt * value
            errors.append(abs(value - math.exp(-1)))
        orders = observed_order(errors)
        self.assertTrue(all(0.95 < p < 1.1 for p in orders))
        with self.assertRaises(ValueError):
            observed_order([0.1, 0])


class StableSegment(unittest.TestCase):
    def test_zero_transfer_all_dissipation_leaves_with_water(self):
        temperature, heat = temperate_segment(10, 0, 2, 5)
        self.assertEqual(temperature, 2.5)
        self.assertEqual(heat, 0)

    def test_strong_transfer_and_energy_closure(self):
        temperature, heat = temperate_segment(10, 1000, 2, 5)
        self.assertAlmostEqual(temperature, .005)
        self.assertAlmostEqual(heat + 10 * temperature, 25)

    def test_weak_transfer_agrees_with_independent_decimal_formula(self):
        from decimal import Decimal, localcontext
        for exposure in (1e-15, 1e-12, 1e-6, .1, 10.):
            with localcontext() as context:
                context.prec = 60
                a = Decimal(str(exposure))
                transfer = 1 - (-a).exp()
                expected = 20 * transfer + 5 * (1 - transfer / a)
            _, heat = temperate_segment(10, 10 * exposure, 2, 5)
            self.assertGreaterEqual(heat, 0)
            self.assertAlmostEqual(heat / float(expected), 1, places=12)

    def test_cold_inlet_tiny_transfer_retains_nonnegative_heat(self):
        temperature, heat = temperate_segment(10, 1e-14, 0, 5)
        self.assertGreater(heat, 0)
        self.assertAlmostEqual(heat / 2.5e-15, 1, places=12)
        self.assertAlmostEqual(temperature, .5)

    def test_invalid_or_refreezing_inputs_rejected(self):
        for values in ((0, 1, 2, 0), (1, -1, 2, 0), (1, 1, -1, 0),
                       (1, 1, 2, -1), (math.inf, 1, 2, 0)):
            with self.assertRaises(ValueError):
                temperate_segment(*values)


class CommandLineContract(unittest.TestCase):
    def test_missing_cases_fail_and_reference_names_are_not_solver_exports(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            generate(root)
            self.assertEqual(len(list(root.glob("*.reference.csv"))), 7)
            self.assertEqual(evaluate(root, root / "report.json", .01, 1e-8), 1)
            result = json.loads((root / "report.json").read_text())
            self.assertFalse(result["passed"])

    def test_cli_identity_fixture_and_corruption(self):
        # These are deliberately copied analytical fixtures to test I/O only.
        # They are NOT results from an independent hydraulic or thermal solver.
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for case in analytical_cases():
                rows = case.rows()
                with (root / (case.name + ".csv")).open("w", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                    writer.writeheader()
                    writer.writerows(rows)
            self.assertEqual(evaluate(root, root / "identity.json", .01, 1e-8), 0)
            path = root / "unequal_outlets.csv"
            rows = read_csv(path)
            rows[10]["outlet_2_m3_s"] = 0
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            self.assertEqual(evaluate(root, root / "corrupt.json", .01, 1e-8), 1)

    def test_csv_duplicate_header_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.csv"
            path.write_text("time_s,time_s\n0,0\n1,1\n")
            with self.assertRaises(InvalidExport):
                read_csv(path)


if __name__ == "__main__":
    unittest.main()
