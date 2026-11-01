import unittest
from run import Clock

class ClockTests(unittest.TestCase):
    def test_constant_rate_offset_and_inverse(self):
        c=Clock([1800,5400,9000],[2.,2.,2.])
        self.assertEqual(c.integral(123,1023),1800.)
        self.assertEqual(c.arrival(123,1800.),1023.)
    def test_zero_flux_pause_and_censoring(self):
        c=Clock([1800,5400,9000],[1.,0.,1.])
        self.assertEqual(c.arrival(3000,1200),7800)
        self.assertIsNone(c.arrival(9000,10000))
    def test_missing_hour_and_extrapolation_refused(self):
        with self.assertRaises(ValueError):Clock([1800,9000],[1.,1.])
        with self.assertRaises(ValueError):Clock([1800,5400],[1.,1.]).integral(-1,1)
if __name__=='__main__':unittest.main()
