import datetime as dt,hashlib,json,tempfile,unittest
from pathlib import Path
from gap_reconstructed_temperature import GapReconstructedTemperatureSeries

BASE=dt.datetime(2019,7,16)
def t(s):return BASE+dt.timedelta(seconds=s)

class GapTests(unittest.TestCase):
    def test_constant_missing_patch_exact_and_accounted_separately(self):
        s=GapReconstructedTemperatureSeries([(t(0),2),(t(30),2),(t(60),float('nan')),(t(90),2),(t(120),2)])
        r=s.integrate(t(0),t(120))
        self.assertEqual(r['observed_seconds'],60)
        self.assertEqual(r['reconstructed_seconds'],60)
        self.assertEqual(r['observed_integral_C_s'],120)
        self.assertEqual(r['reconstructed_integral_C_s'],120)
        self.assertEqual(r['integral_C_s'],240)
        self.assertEqual(r['coverage_fraction'],1)

    def test_linear_partial_gap_integral_and_persistence_sensitivities(self):
        s=GapReconstructedTemperatureSeries([(t(0),0),(t(30),3),(t(60),float('nan')),(t(90),9),(t(120),12)])
        r=s.integrate(t(45),t(75))
        self.assertEqual(r['observed_seconds'],0)
        self.assertEqual(r['reconstructed_integral_C_s'],180)
        self.assertEqual(s.integrate(t(45),t(75),method='constant_previous_endpoint')['integral_C_s'],90)
        self.assertEqual(s.integrate(t(45),t(75),method='constant_next_endpoint')['integral_C_s'],270)

    def test_over_24_hours_not_filled_and_no_extrapolation(self):
        s=GapReconstructedTemperatureSeries([(t(0),1),(t(30),1),(t(90000),1),(t(90030),1)])
        r=s.integrate(t(-10),t(90040))
        self.assertEqual(r['observed_seconds'],60)
        self.assertEqual(r['reconstructed_seconds'],0)
        self.assertEqual(r['max_uncovered_gap_seconds'],89970)
        with self.assertRaises(ValueError):GapReconstructedTemperatureSeries([(t(0),1),(t(30),1)],max_reconstruction_gap_seconds=86401)

    def test_raw_final_offset_applied_once_for_observed_and_reconstructed(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'products').mkdir()
            p=root/'products/temp3_corr_raw.csv';p.write_text('Date_Time,Temperature\n2019-07-16T0:0:0,0.1\n2019-07-16T0:0:30,0.1\n2019-07-16T0:1:0,NaN\n2019-07-16T0:1:30,0.1\n')
            (root/'SOURCE_MANIFEST.json').write_text(json.dumps({'members':[{'path':'products/temp3_corr_raw.csv','verified_crc32':True,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}]}))
            s=GapReconstructedTemperatureSeries.from_archive(root,3)
            r=s.integrate(t(0),t(90))
            self.assertAlmostEqual(r['integral_C_s'],.135*90)
            self.assertAlmostEqual(r['observed_integral_C_s'],.135*30)
            self.assertAlmostEqual(r['reconstructed_integral_C_s'],.135*60)

    def test_positive_clamp_exact_within_reconstructed_linear_gap(self):
        s=GapReconstructedTemperatureSeries([(t(0),-2),(t(30),-1),(t(60),float('nan')),(t(90),1),(t(120),2)])
        r=s.integrate_positive(t(30),t(90))
        self.assertEqual(r['integral_C_s'],15)
        self.assertEqual(r['unclipped_signed_integral_C_s'],0)

if __name__=='__main__':unittest.main()
