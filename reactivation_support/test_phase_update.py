import unittest
import numpy as np
from .phase_update import update

class PhaseTests(unittest.TestCase):
    def test_closed_box_total_mass_and_energy(self):
        h=.2;T=2.;H=1.;dz=.001;ri,rw,cp,L=917.,1000.,4217.,334000.
        r=update(0.,h,0.,0.,h*T,dz,liquid_energy_delta_J_m2=-ri*L*dz)
        self.assertAlmostEqual(rw*r['h']+ri*(H-dz),rw*h+ri*H,places=10)
        self.assertAlmostEqual(rw*cp*r['hT']+ri*L*dz,rw*cp*h*T,places=8)
        self.assertAlmostEqual(float(r['required_external_energy_J_m2']),0.,places=12)
    def test_equal_densities_leave_surface_stationary(self):
        r=update(-1.,.2,0.,0.,.1,.03,liquid_energy_delta_J_m2=0.,rho_i=1000.,rho_w=1000.)
        self.assertAlmostEqual(float(r['zb']+r['h']),-.8,places=12)
    def test_ice_water_density_jump_surface_falls_less_than_floor(self):
        r=update(0.,.2,0.,0.,.1,.01,liquid_energy_delta_J_m2=0.)
        self.assertAlmostEqual(float(r['zb']+r['h']-.2),-.00083,places=12)
    def test_added_cold_water_dilutes_without_second_heat_sink(self):
        r=update(0.,.2,.1,-.02,.1,.01,liquid_energy_delta_J_m2=0.)
        self.assertEqual(float(r['hT']),.1)
        self.assertLess(float(r['temperature']),.5)
        self.assertEqual(float(r['hu']),.1)
        self.assertEqual(float(r['hv']),-.02)
    def test_source_momentum_explicit(self):
        r=update(0.,.2,0.,0.,.1,.01,liquid_energy_delta_J_m2=0.,melt_u=2.)
        self.assertAlmostEqual(float(r['hu']),.01834,places=12)
    def test_no_melt_identity(self):
        r=update(0.,.2,.1,0.,.1,0.,liquid_energy_delta_J_m2=0.)
        self.assertEqual(float(r['h']),.2);self.assertEqual(float(r['zb']),0.)
    def test_energy_overdraw_rejected(self):
        with self.assertRaises(ValueError):update(0.,.1,0.,0.,0.,.1,liquid_energy_delta_J_m2=-1.)
    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):update(0.,.1,0.,0.,.1,float('nan'),liquid_energy_delta_J_m2=0.)
    def test_array_input_not_mutated(self):
        h=np.full((3,4),.2);old=h.copy()
        update(np.zeros_like(h),h,h*0,h*0,h*.5,.001,liquid_energy_delta_J_m2=0.)
        np.testing.assert_array_equal(h,old)
    def test_randomized_mass_and_sensible_balance(self):
        rng=np.random.default_rng(1101)
        h=rng.uniform(.01,1.,(9,7));T=rng.uniform(.1,2.,h.shape);dz=rng.uniform(0.,.001,h.shape)
        r=update(h*0,h,h*.1,h*0,h*T,dz,liquid_energy_delta_J_m2=0.)
        np.testing.assert_allclose(1000*(r['h']-h),917*dz,rtol=1e-11,atol=2e-13)
        np.testing.assert_allclose(1000*4217*(r['hT']-h*T),0.,atol=1e-10)

if __name__=='__main__':unittest.main()
