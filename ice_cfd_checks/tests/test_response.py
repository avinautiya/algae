import math
import unittest

from ice_cfd_checks.linear_response import analyze, exponential, largest_symmetric_eigenvalue, response_gain


class ResponseTests(unittest.TestCase):
    def test_diagonal_decay(self):
        self.assertAlmostEqual(response_gain([[-1,0],[0,-2]],0.7),math.exp(-1.4),places=13)

    def test_jordan_block_against_closed_form(self):
        time=0.8
        m=exponential([[-1,3],[0,-1]],time)
        self.assertAlmostEqual(m[0][0],math.exp(-time),places=13)
        self.assertAlmostEqual(m[0][1],3*time*math.exp(-time),places=13)
        self.assertAlmostEqual(m[1][1],math.exp(-time),places=13)

    def test_non_normal_exponential_and_gain_against_two_by_two_formula(self):
        t=math.log(2); a,d=math.exp(-t),math.exp(-2*t); b=10*(a-d)
        m=exponential([[-1,10],[0,-2]],t)
        self.assertAlmostEqual(m[0][1],b,places=13)
        trace=a*a+b*b+d*d; determinant=(a*d)**2
        exact=(trace+math.sqrt(trace*trace-4*determinant))/2
        self.assertAlmostEqual(response_gain([[-1,10],[0,-2]],t),exact,places=12)
        self.assertGreater(exact,1)

    def test_orthogonal_rotation_preserves_norm(self):
        self.assertAlmostEqual(response_gain([[0,-1],[1,0]],1.8),1,places=13)

    def test_input_output_maps_and_jacobi(self):
        self.assertAlmostEqual(response_gain([[-1,10],[0,-2]],math.log(2),[[0],[1]],[[1,0]]),2.5**2,places=12)
        self.assertAlmostEqual(largest_symmetric_eigenvalue([[2,1,0],[1,2,0],[0,0,-5]]),3,places=12)

    def test_normalization_and_provenance_required(self):
        with self.assertRaises(ValueError):
            analyze({"A_per_s":[[-1]],"times_s":[0,1]})
        with self.assertRaises(ValueError):
            exponential([[-1]],-1)


if __name__=="__main__":
    unittest.main()
