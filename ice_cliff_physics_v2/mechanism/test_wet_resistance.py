import unittest,math
from wet_resistance import prediction,thickness_limit
class ResistanceTests(unittest.TestCase):
    def test_series_network(self):
        r=prediction(water_capacity_W_K=1e12,film_conductance_W_K=100,ice_conductance_W_K=20,inlet_superheat_K=.5)
        self.assertAlmostEqual(r['covered_ice_power_W'],.5/(1/100+1/20),places=8)
    def test_finite_heat_never_exceeded(self):
        for F in [.1,10,1000]:
            r=prediction(water_capacity_W_K=F,film_conductance_W_K=100,ice_conductance_W_K=20,inlet_superheat_K=.5)
            self.assertLessEqual(r['covered_ice_power_W'],r['bare_ice_power_W'])
            self.assertLessEqual(r['bare_ice_power_W'],F*.5)
    def test_thickness_boundary_exact(self):
        A=.2;k=1;F=200;G=100;target=.8
        d=thickness_limit(wet_area_m2=A,conductivity_W_m_K=k,water_capacity_W_K=F,film_conductance_W_K=G,required_fraction=target)
        r=prediction(water_capacity_W_K=F,film_conductance_W_K=G,ice_conductance_W_K=k*A/d,inlet_superheat_K=.5)
        self.assertAlmostEqual(r['thermal_reconnection_fraction'],target,places=12)
    def test_unit_recovery_requires_zero_thickness(self):
        self.assertEqual(thickness_limit(wet_area_m2=.2,conductivity_W_m_K=1,water_capacity_W_K=200,film_conductance_W_K=100,required_fraction=1),0)
    def test_rewetting_can_cool_without_cold_water(self):
        r=prediction(water_capacity_W_K=200,film_conductance_W_K=100,ice_conductance_W_K=20,inlet_superheat_K=.5,absorbed_external_heat_W=100)
        self.assertLess(r['wet_minus_dry_ice_power_W'],0)
        self.assertLess(r['water_to_debris_heat_W'],0)
        self.assertGreater(r['outgoing_water_temperature_C'],.5)
        self.assertLess(abs(r['energy_residual_W']),1e-12)
    def test_rewetting_can_warm(self):
        r=prediction(water_capacity_W_K=200,film_conductance_W_K=100,ice_conductance_W_K=20,inlet_superheat_K=.5,absorbed_external_heat_W=5)
        self.assertGreater(r['wet_minus_dry_ice_power_W'],0)
        self.assertGreater(r['water_to_debris_heat_W'],0)
    def test_crossover_exact(self):
        for G in [.1,100,100000]:
            r=prediction(water_capacity_W_K=200,film_conductance_W_K=G,ice_conductance_W_K=20,inlet_superheat_K=.5,absorbed_external_heat_W=10)
            self.assertAlmostEqual(r['wet_minus_dry_ice_power_W'],0,places=11)
            self.assertAlmostEqual(r['water_to_debris_heat_W'],0,places=11)
    def test_invalid(self):
        with self.assertRaises(ValueError):prediction(water_capacity_W_K=0,film_conductance_W_K=100,ice_conductance_W_K=20,inlet_superheat_K=.5)
        with self.assertRaises(ValueError):thickness_limit(wet_area_m2=.2,conductivity_W_m_K=1,water_capacity_W_K=200,film_conductance_W_K=100,required_fraction=1.1)
if __name__=='__main__':unittest.main()
