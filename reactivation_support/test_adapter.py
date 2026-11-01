import unittest
import numpy as np
from .snapshot.swe2d import Sim
from .integration_adapter import conservative_lower_bed

class AdapterTests(unittest.TestCase):
    def test_phase_source_on_pinned_snapshot(self):
        s=Sim(dx=.1,Q=0.,closed_out=True);s.last_qacc=np.where(s.ice,100.,0.)
        before=float(s.h.sum()*s.dx**2)
        v,_,_=conservative_lower_bed(s,3600.)
        added=float(s.h.sum()*s.dx**2)-before
        self.assertAlmostEqual(added,v*917/1000,places=12)
        self.assertAlmostEqual(s.phase_water_added_m3,added,places=12)
    def test_nondefault_materials_consistent(self):
        s=Sim(dx=.1,Q=0.,closed_out=True);s.last_qacc=np.where(s.ice,100.,0.)
        before=float(s.h.sum()*s.dx**2)
        v,_,_=conservative_lower_bed(s,3600.,rho_i=900.,rho_w=1100.,latent=300000.)
        self.assertAlmostEqual(float(s.h.sum()*s.dx**2)-before,v*900/1100,places=12)
if __name__=='__main__':unittest.main()
