import datetime as dt
import unittest
from native_temperature import TemperatureSeries

BASE=dt.datetime(2019,7,16)
def t(s):return BASE+dt.timedelta(seconds=s)

class NativeIntegralTests(unittest.TestCase):
    def test_exact_linear_integral_at_arbitrary_partial_segment_endpoints(self):
        s=TemperatureSeries([(t(0),1.),(t(10),3.),(t(30),7.)],60)
        result=s.integrate(t(5),t(25))
        self.assertAlmostEqual(result['integral_C_s'],80.)
        self.assertEqual(result['coverage_fraction'],1.)
        self.assertEqual(result['max_uncovered_gap_seconds'],0.)

    def test_gap_and_nan_are_not_bridged(self):
        s=TemperatureSeries([(t(0),1.),(t(10),1.),(t(20),float('nan')),(t(30),2.),(t(40),2.),(t(200),3.),(t(210),3.)],60)
        r=s.integrate(t(0),t(210))
        self.assertEqual(r['covered_seconds'],30.)
        self.assertEqual(r['integral_C_s'],60.)
        self.assertEqual(r['max_uncovered_gap_seconds'],160.)

    def test_irregular_samples_weight_time_not_sample_count(self):
        s=TemperatureSeries([(t(0),0.),(t(1),1.),(t(2),2.),(t(10),10.)],60)
        self.assertEqual(s.integrate(t(0),t(10))['integral_C_s'],50.)
        self.assertNotEqual(50./10,sum([0,1,2,10])/4)

    def test_sensitivity_offset_only_on_supported_time(self):
        s=TemperatureSeries([(t(0),1.),(t(10),1.)],60)
        r=s.integrate(t(-5),t(15),extra_offset_C=.05)
        self.assertEqual(r['covered_seconds'],10.)
        self.assertAlmostEqual(r['integral_C_s'],10.5)
        self.assertEqual(r['max_uncovered_gap_seconds'],5.)

    def test_duplicate_or_backward_time_refused(self):
        with self.assertRaises(ValueError):TemperatureSeries([(t(0),1.),(t(0),1.)])
        with self.assertRaises(ValueError):TemperatureSeries([(t(1),1.),(t(0),1.)])

    def test_timezone_must_be_reconciled_upstream(self):
        s=TemperatureSeries([(t(0),1.),(t(10),1.)])
        with self.assertRaises(ValueError):s.integrate(t(0).replace(tzinfo=dt.timezone.utc),t(10))

    def test_positive_integral_resolves_zero_crossing_exactly(self):
        s=TemperatureSeries([(t(0),-1.),(t(10),1.)])
        r=s.integrate_positive(t(0),t(10))
        self.assertEqual(r['integral_C_s'],2.5)
        self.assertEqual(r['unclipped_signed_integral_C_s'],0.)
        self.assertEqual(s.integrate_positive(t(2.5),t(7.5))['integral_C_s'],.625)
        self.assertEqual(s.integrate_positive(t(0),t(10),extra_offset_C=-2)['integral_C_s'],0.)

if __name__=='__main__':unittest.main()
