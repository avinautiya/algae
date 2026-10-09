import unittest
from agency_design.prombio import observation, load_2024


class PrombioTests(unittest.TestCase):
    def test_censoring_not_zero(self):
        r = observation("<1000")
        self.assertEqual(r["kind"], "left_censored")
        self.assertIsNone(r["value"])
        self.assertEqual(r["upper_bound"], 1000)

    def test_decimal_comma_and_zero(self):
        self.assertEqual(observation("4,00E+04")["value"], 40000)
        self.assertEqual(observation("0")["kind"], "observed")
        self.assertEqual(observation("NA")["kind"], "missing")
        self.assertIsNone(observation("no counts")["value"])

    def test_invalid(self):
        for value in ("-1", "inf", "<0"):
            with self.assertRaises(ValueError):
                observation(value)

    def test_real_table_checksum_and_inventory(self):
        records, provenance = load_2024()
        self.assertEqual(len(records), 118)
        self.assertEqual(sum(r["ice_cells_ml"]["kind"] == "observed" for r in records), 108)
        self.assertEqual(len({r["station_area"] for r in records}), 5)
        self.assertEqual(sum(r["sampling_stratum"] == "under radiometer" for r in records), 31)
        self.assertEqual(sum(r["sampling_stratum"] == "randomized grid" for r in records), 75)
        self.assertTrue(all(r["latitude"] is None for r in records))
        self.assertEqual(provenance["license"], "CC0 1.0")


if __name__ == "__main__":
    unittest.main()
