import math
import unittest
from uncertainty import estimate


class Tests(unittest.TestCase):
    def test_exact_known_power_law(self):
        hs=[.1,.05,.025];p=.8;truth=2.18;vals=[truth+3*h**p for h in hs]
        r=estimate(vals,hs)
        self.assertAlmostEqual(r['observed_order'],p)
        self.assertAlmostEqual(r['extrapolated_value'],truth)

    def test_published_rounded_values(self):
        r=estimate([3.165,2.721,2.478],[.1,.05,.025])
        self.assertAlmostEqual(r['observed_order'],math.log(.444/.243)/math.log(2))
        self.assertAlmostEqual(r['gci_absolute'],.3672201492537326)
        self.assertLess(r['extrapolated_value'],2.19)

    def test_sign_and_datum_invariance(self):
        a=estimate([3.,2.,1.5],[.1,.05,.025])
        b=estimate([-3.,-2.,-1.5],[.1,.05,.025])
        self.assertAlmostEqual(a['observed_order'],b['observed_order'])
        self.assertAlmostEqual(a['gci_absolute'],b['gci_absolute'])

    def test_invalid_or_censored(self):
        for values in ([1.,None,2.],[1.,math.inf,2.],[True,1.,2.]):
            with self.assertRaises(ValueError):estimate(values,[.1,.05,.025])
        with self.assertRaises(ValueError):estimate([3.,2.,1.],[.1,.05,.02])

    def test_refuses_false_precision(self):
        for vals,status in [([1.,1.,1.],'zero_difference'),([3.,1.,2.],'oscillatory'),([4.,3.,1.],'noncontracting')]:
            r=estimate(vals,[.1,.05,.025]);self.assertIn(status,r['status'])
            self.assertNotIn('gci_absolute',r)


if __name__=='__main__':unittest.main()
