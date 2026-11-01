import unittest
from audit_head_drop import head_drop, f_nominal, f_stage_fixed

class HeadDropTests(unittest.TestCase):
    def test_equal_depths_cancel(self):
        self.assertAlmostEqual(head_drop(5, .6, 4, .6), 1)

    def test_equal_heads_give_zero(self):
        self.assertAlmostEqual(head_drop(4, .6, 4, .6), 0)
        self.assertAlmostEqual(head_drop(5, .1, 4, 1.1), 0)

    def test_common_elevation_datum_cancels(self):
        self.assertAlmostEqual(head_drop(5, .6, 4, .4), head_drop(1005, .6, 1004, .4))

    def test_exchange_endpoints_reverses_sign(self):
        self.assertAlmostEqual(head_drop(5,.6,4,.4), -head_drop(4,.4,5,.6))

    def test_negative_depth_refused(self):
        with self.assertRaises(ValueError):head_drop(5,-.6,4,.4)

    def test_friction_scales_with_drop_in_both_branches(self):
        self.assertAlmostEqual(f_nominal(1,253.7,1,.8)/f_nominal(2,253.7,1,.8), .5)
        self.assertAlmostEqual(f_stage_fixed(1,253.7,1,.8,.6,.4)/f_stage_fixed(2,253.7,1,.8,.6,.4), .5)

if __name__=='__main__':unittest.main()
