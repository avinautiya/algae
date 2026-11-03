import math, unittest
from energy_contract import *
class Checks(unittest.TestCase):
    def test_radiometer(self):
        eps=.96; incoming=300.; T=290.
        outgoing=eps*SIGMA*T**4+(1-eps)*incoming
        self.assertAlmostEqual(radiometric_temperature_C(outgoing,incoming,eps,surface_is_bare_debris=True),T-273.15,11)
    def test_phase_refusal(self):
        for bare in [False,None]:
            with self.assertRaises(ValueError): radiometric_temperature_C(350,300,.96,surface_is_bare_debris=bare)
    def test_steady_linear_profile(self):
        # Independent Fourier gradient gives q=k*Ts/d. No storage change.
        b=conduction_energy_bounds(.96,.06,1.8e6,10*3600,0,10)
        self.assertAlmostEqual(b['quasisteady_signed_energy_J_m2'],.96*10/.06*3600)
        self.assertLessEqual(b['signed_energy_lower_J_m2'],b['quasisteady_signed_energy_J_m2'])
    def test_weighted_mode_identity(self):
        # T=A sin(pi*z/d) exp(-alpha*pi²*t/d²), surface zero.
        # I=A*d/pi*exp(...); -C I_t equals independent basal Fourier flux.
        d=.06;k=.96;C=1.8e6;A=5.;t=100.;alpha=k/C
        Iprime=-A*d/math.pi*(alpha*math.pi**2/d**2)*math.exp(-alpha*math.pi**2*t/d**2)
        flux=k*A*math.pi/d*math.exp(-alpha*math.pi**2*t/d**2)
        self.assertAlmostEqual(-C*Iprime,flux,10)
    def test_storage_bound(self):
        b=conduction_energy_bounds(.96,.06,1.8e6,0,0,30)
        self.assertEqual(b['maximum_storage_correction_J_m2'],1620000)
        self.assertLess(b['signed_energy_lower_J_m2'],0)
    def test_melt_scope_and_units(self):
        with self.assertRaises(ValueError): signed_energy_to_net_melt_mmwe(334000,temperate_ice_confirmed=False,all_other_interface_terms_zero_confirmed=True)
        self.assertEqual(signed_energy_to_net_melt_mmwe(334000,temperate_ice_confirmed=True,all_other_interface_terms_zero_confirmed=True),1)
    def test_basal_temperature_in_range(self):
        with self.assertRaises(ValueError): conduction_energy_bounds(.96,.06,1.8e6,1,5,30)
    def test_invalid(self):
        for k in [0,-1,float('nan')]:
            with self.assertRaises(ValueError): conduction_energy_bounds(k,.06,1,1,0,30)
if __name__=='__main__': unittest.main()
