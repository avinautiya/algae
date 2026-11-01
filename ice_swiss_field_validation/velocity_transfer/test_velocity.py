import math
import unittest
from run_velocity_test import diameter, velocity_basis, fit_amplitude, metrics

class VelocityChecks(unittest.TestCase):
    def test_diameter_from_independent_geometry_only(self):
        self.assertAlmostEqual(diameter(2,.6,.4),4/3)

    def test_no_target_or_discharge_in_diameter(self):
        import inspect
        self.assertEqual(list(inspect.signature(diameter).parameters),['width','up_depth','down_depth'])

    def test_fit_recovers_known_coefficient(self):
        amplitude=fit_amplitude([1.,2.,3.],[2.,4.,6.])
        self.assertAlmostEqual(amplitude,2.)
        self.assertAlmostEqual(1/amplitude**2,.25)

    def test_head_scaling(self):
        self.assertAlmostEqual(velocity_basis(4,1,100)/velocity_basis(1,1,100),2)

    def test_invalid_inputs_refused(self):
        for head in [0,-1,float('nan')]:
            with self.assertRaises(ValueError):velocity_basis(head,1,100)

    def test_metric_arithmetic(self):
        m=metrics([1.,2.],[2.,0.])
        self.assertEqual(m['MAE_m_s'],1.5)
        self.assertAlmostEqual(m['RMSE_m_s'],math.sqrt(2.5))

if __name__=='__main__':unittest.main()
