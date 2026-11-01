import datetime as dt
import unittest
from run_audit import block_errors,TemperatureSeries

BASE=dt.datetime(2019,7,16)
def t(s):return BASE+dt.timedelta(seconds=s)

class MaskTests(unittest.TestCase):
    def test_constant_signal_all_methods_exact(self):
        s=TemperatureSeries([(t(0),.2),(t(30),.2),(t(60),.2)])
        rows=block_errors(s,t(0),t(60))
        for r in rows:self.assertAlmostEqual(r['signed_integral_error_C_s'],0)

    def test_linear_signal_bridge_exact_constant_endpoints_bracket(self):
        s=TemperatureSeries([(t(0),0),(t(30),3),(t(60),6)])
        rows={r['method']:r for r in block_errors(s,t(0),t(60))}
        self.assertEqual(rows['linear_endpoint_bridge']['signed_integral_error_C_s'],0)
        self.assertEqual(rows['constant_previous_endpoint']['signed_integral_error_C_s'],-180)
        self.assertEqual(rows['constant_next_endpoint']['signed_integral_error_C_s'],180)

    def test_partially_supported_block_excluded(self):
        s=TemperatureSeries([(t(0),0),(t(30),1),(t(500),2),(t(530),3)])
        self.assertIsNone(block_errors(s,t(0),t(530)))

if __name__=='__main__':unittest.main()
