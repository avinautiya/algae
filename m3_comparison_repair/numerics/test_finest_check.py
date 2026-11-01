import types
import unittest
import numpy as np
from finest_check import operator_preflight,plan


class Fake:
    def __init__(self,correct):
        self.zb=np.array([[-2.,-.8],[0.,0.]])
        self.ice=np.array([[False,True],[True,False]])
        self.pool_floor=-1.;self.dx=.1;self.correct=correct
    def lower_bed(self,dt):
        d=self.last_qacc*dt/(917*334000)*self.ice
        if self.correct:
            before=self.zb.copy()
            self.zb=np.where(self.ice,np.maximum(self.zb-d,self.pool_floor),self.zb)
            vol=float(((before-self.zb)*self.ice).sum()*self.dx**2)
        else:
            self.zb=np.maximum(self.zb-d,self.pool_floor)
            vol=float(d.sum()*self.dx**2)
        return vol,vol,0.


class Tests(unittest.TestCase):
    def test_bad_operator_blocks(self):
        c=types.SimpleNamespace(build=lambda *a,**k:Fake(False),RHO_I=917,LF=334000)
        r=operator_preflight(c)
        self.assertFalse(r['passed']);self.assertEqual(r['zero_dt']['nonice_changed_cells'],1)
        self.assertGreater(r['cap_activation']['reported_minus_actual_m3'],0.)

    def test_correct_operator_passes_local_checks_only(self):
        c=types.SimpleNamespace(build=lambda *a,**k:Fake(True),RHO_I=917,LF=334000)
        self.assertTrue(operator_preflight(c)['passed'])

    def test_cost_scope(self):
        p=plan(3.5,None)
        self.assertAlmostEqual(p['estimated_compute_wall_s'],2704.)
        self.assertIn('early contact-loss',p['supports'])
        self.assertFalse(p['checkpoint_support'])


if __name__=='__main__':unittest.main()
