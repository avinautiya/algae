import itertools,unittest
from measurement_design import reading_resolution

class DesignTests(unittest.TestCase):
    def test_threshold_has_adversarial_tie(self):
        fa,fb,x,y=.27,.43,846.,315.
        d=reading_resolution(fa,fb,x,y)['critical_bound_mmwe']
        vals=[]
        for es in itertools.product((-d,d),repeat=3):
            xx=x+es[1]-es[0];yy=y+es[2]-es[1]
            vals.append((fa*xx-yy)**2-(fb*xx-yy)**2)
        self.assertAlmostEqual(min(abs(v) for v in vals),0,places=8)
    def test_invalid(self):
        with self.assertRaises(ValueError):reading_resolution(-1,1,1,1)
        self.assertEqual(reading_resolution(1,1,2,1)['status'],'identical_predictions')
