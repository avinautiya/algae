import unittest,math
from model import Case
from rescue_certificate import *

class CertificateTests(unittest.TestCase):
    def test_partial_submergence_force_balance(self):
        c=Case();D=.2;h=cube_critical_depth(D,c.friction,c.drag,c.grain_density/c.water_density)
        U=math.sqrt(c.gravity*h)
        drag=.5*c.water_density*c.drag*D*h*U**2
        friction=c.friction*c.gravity*D**2*(c.grain_density*D-c.water_density*h)
        self.assertAlmostEqual(drag/friction,1.,places=12)
        self.assertLess(h,D)
    def test_full_submergence_force_balance(self):
        c=Case(friction=1.);D=.2
        h=cube_critical_depth(D,c.friction,c.drag,c.grain_density/c.water_density)
        self.assertGreater(h,D)
        drag=.5*c.water_density*c.drag*D**2*c.gravity*h
        friction=c.friction*c.gravity*D**3*(c.grain_density-c.water_density)
        self.assertAlmostEqual(drag/friction,1.,places=12)
    def test_zero_friction_threshold_reduces_to_rewetting(self):
        c=Case(friction=0.);g=.2
        K=(2/3)**1.5*math.sqrt(c.gravity)*c.width
        self.assertAlmostEqual(minimum_sliding_discharge(c,g,.2),K*g**1.5,places=12)
    def test_longer_flood_weakens_exclusion_conservatively(self):
        c=Case()
        values=[]
        for tau in (3600.,21600.,86400.):
            r=exclude_rescue(c,toe_ice_z=-.1,bypass_ice_z=-.5,
                             protected_solid_stock_m3=.024,grain_D_m=.2,
                             peak_discharge_m3_s=.14,duration_s=tau)
            values.append(r['optimistic_minimum_sliding_discharge_m3_s'])
        self.assertGreater(values[0],values[1]);self.assertGreater(values[1],values[2])
    def test_bypass_incision_raises_rescue_threshold(self):
        c=Case()
        self.assertGreater(minimum_sliding_discharge(c,.3,.2),minimum_sliding_discharge(c,.1,.2))
if __name__=='__main__':unittest.main()

class ThermalLockoutTests(unittest.TestCase):
    def test_extra_heat_can_only_weaken_exclusion(self):
        c=Case()
        a=thermal_lockout(c,toe_ice_z=-.1,bypass_ice_z=-.5,peak_discharge_m3_s=.2,duration_s=21600)
        b=thermal_lockout(c,toe_ice_z=-.1,bypass_ice_z=-.5,peak_discharge_m3_s=.2,duration_s=21600,extra_toe_heat_W_m2=100.)
        self.assertGreater(a['optimistic_remaining_ice_gap_m'],b['optimistic_remaining_ice_gap_m'])
    def test_zero_duration_recovers_exact_empty_route_threshold(self):
        c=Case();K=(2/3)**1.5*math.sqrt(c.gravity)*c.width
        head=(c.discharge/K)**(2/3)
        for factor in (.99,1.01):
            r=thermal_lockout(c,toe_ice_z=factor*head,bypass_ice_z=0,peak_discharge_m3_s=.2,duration_s=0)
            self.assertEqual(r['reconnection_excluded_even_after_perfect_clearing'],factor>1)
