import math
import unittest
from dataclasses import replace
import numpy as np
from model import Case,evolve,partition

class RoutingTests(unittest.TestCase):
    def test_symmetric_split_and_discharge_ledger(self):
        q1,q2,eta=partition(np.array([.1,.2]),np.zeros(2),np.zeros(2),np.ones(2),np.ones(2))
        np.testing.assert_allclose(q1,[.05,.1],rtol=1e-12)
        np.testing.assert_allclose(q1+q2,[.1,.2],rtol=1e-12)
    def test_dry_threshold_is_geometric(self):
        Q,K=.1,1.
        for multiplier in (.99,1.01):
            q1,q2,eta=partition(Q,multiplier*(Q/K)**(2/3),0,K,K)
            self.assertEqual(q1>0,multiplier<1)
    def test_no_export_no_heat_exact_stock_and_abandonment(self):
        c=Case(h_water=0.,allow_export=False)
        r=evolve([c],days=10,dt=600)
        As=(1-c.porosity)*c.width*c.length
        S0=c.initial_debris_depth*As
        source=(1-c.porosity)*c.mantle_thickness*c.source_width*c.face_heat_W_m2/(c.ice_density*c.latent)
        target=(c.initial_bypass_z+(.1/((2/3)**1.5*math.sqrt(c.gravity)*c.width))**(2/3)-c.initial_toe_ice_z)*As
        exact=(target-S0)/source/86400
        self.assertLessEqual(abs(r['first_dry_day'][0]-exact),600/86400+1e-10)
        self.assertAlmostEqual(sum(r['final_stock_m3'][0]),S0+source*10*86400,places=10)
    def test_energy_and_solid_ledgers(self):
        r=evolve([Case(),Case(source_width=0),Case(allow_ice_incision=False)],days=2,dt=600)
        self.assertLess(r['max_solid_ledger_residual_m3'],1e-12)
        np.testing.assert_allclose(np.array(r['energy_input_J'])-r['energy_advected_out_J'],r['latent_J'],rtol=2e-12,atol=1e-6)
    def test_no_heat_and_no_supply_has_no_synthetic_ice_loss(self):
        r=evolve([Case(h_water=0,allow_supply=False)],days=2,dt=600)
        self.assertEqual(r['toe_melt_m3'][0],0)
        self.assertEqual(r['bypass_melt_m3'][0],0)
    def test_debris_stocks_nonnegative(self):
        r=evolve([Case(friction=0.,length=.01)],days=.1,dt=600)
        self.assertGreaterEqual(min(r['final_stock_m3'][0]),0)
        self.assertLess(r['max_solid_ledger_residual_m3'],1e-12)
    def test_invalid_material(self):
        for kwargs in ({'porosity':1.},{'discharge':0.},{'water_temperature':-1.},{'source_fractions':(.2,.2,.2)}):
            with self.assertRaises(ValueError):Case(**kwargs)
    def test_elevation_includes_melt_accommodation(self):
        c=Case(allow_supply=False,allow_export=False)
        r=evolve([c],days=1,dt=300)
        expected=c.initial_toe_ice_z+c.initial_debris_depth-r['toe_melt_m3'][0]/(c.width*c.length)
        self.assertAlmostEqual(r['final_toe_z'][0],expected,places=11)

if __name__=='__main__':unittest.main()

class CoverageTests(unittest.TestCase):
    def test_monolayer_parallel_heat_flux(self):
        c=Case(source_fractions=(0.,0.,1.),allow_supply=False,allow_export=False,footprint_shield=True)
        duration=1.
        r=evolve([c],days=duration/86400,dt=duration)
        K=(2/3)**1.5*math.sqrt(c.gravity)*c.width
        q,_,_=partition(c.discharge,c.initial_toe_ice_z+c.initial_debris_depth,c.initial_bypass_z,K,K)
        S=c.initial_debris_depth*(1-c.porosity)*c.width*c.length
        cover=min(S/(.2*c.width*c.length),1.)
        h=(1-cover)*c.h_water+cover/(1/c.h_water+.2/c.debris_k)
        cpq=c.water_density*c.water_cp*q
        P=cpq*c.water_temperature*(-math.expm1(-h*c.width*c.length/cpq))
        expected=P*duration/(c.ice_density*c.latent)
        self.assertAlmostEqual(r['toe_melt_m3'][0]/expected,1.,places=11)
    def test_patchy_and_uniform_collapse_to_bare_without_grains(self):
        a=Case(initial_debris_depth=0,allow_supply=False)
        r=evolve([a,replace(a,footprint_shield=False)],days=.01,dt=60)
        self.assertAlmostEqual(r['toe_melt_m3'][0],r['toe_melt_m3'][1],places=12)
    def test_extra_heat_is_in_energy_ledger(self):
        r=evolve([Case(extra_toe_heat_W_m2=100.)],days=.1,dt=60)
        self.assertGreater(r['external_toe_heat_J'][0],0)
        residual=r['energy_input_J'][0]-r['energy_advected_out_J'][0]-r['latent_J'][0]
        self.assertLess(abs(residual)/r['energy_input_J'][0],1e-12)
