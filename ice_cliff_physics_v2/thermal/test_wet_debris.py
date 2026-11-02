import math
import unittest
from wet_debris import ThermalState as S, ThermalParameters as P, UnsupportedThermalState, advance


class WetDebrisTests(unittest.TestCase):
    def close_ledger(self,r):
        scale=max(1.,abs(r.water_inlet_heat_J),abs(r.atmospheric_heat_J),abs(r.debris_storage_change_J))
        self.assertLess(abs(r.global_energy_residual_J)/scale,2e-12)
        self.assertLess(abs(r.water_exchange_residual_J)/scale,2e-12)
        self.assertGreaterEqual(r.passivity_margin_J,-2e-12*scale)

    def test_bare_only_finite_heat(self):
        r=advance(S(0),P(0,0,200,0),30,.2,100)
        expected=100*.2*(-math.expm1(-2))*30
        self.assertAlmostEqual(r.bare_ice_heat_J,expected)
        self.assertAlmostEqual(r.outlet_mean_temperature_C,.2*math.exp(-2))
        self.assertEqual(r.covered_ice_heat_J,0)
        self.close_ledger(r)

    def test_dry_conduction_and_atmosphere(self):
        r=advance(S(4),P(2000,10,300,400),100,0,0,20)
        expected=2+(4-2)*math.exp(-.5)
        self.assertAlmostEqual(r.state.temperature_C,expected)
        self.assertEqual(r.water_to_debris_heat_J,0)
        self.assertEqual(r.water_to_bare_heat_J,0)
        self.assertIsNone(r.outlet_mean_temperature_C)
        self.close_ledger(r)

    def test_zero_capacity_algebraic_series_resistance(self):
        r=advance(S(5),P(0,20,0,80),15,.5,200)
        H=200*(-math.expm1(-80/200))
        t=H*.5/(H+20)
        self.assertAlmostEqual(r.state.temperature_C,t)
        self.assertAlmostEqual(r.covered_ice_heat_J,20*t*15)
        self.assertEqual(r.debris_storage_change_J,0)
        self.close_ledger(r)

    def test_hot_debris_can_heat_the_water(self):
        r=advance(S(20),P(10000,10,0,100),50,.1,200)
        self.assertLess(r.water_to_debris_heat_J,0)
        self.assertGreater(r.water_outlet_heat_J,r.water_inlet_heat_J)
        self.assertGreater(r.covered_ice_heat_J,0)
        self.close_ledger(r)

    def test_parallel_hot_debris_can_deliver_heat_to_bare_ice(self):
        r=advance(S(10),P(10000,1,100,100),10,0,1e-9)
        self.assertGreater(r.bare_ice_heat_J,0)
        self.assertLess(r.water_to_debris_heat_J,0)
        self.assertLess(r.total_ice_heat_J,100000)
        self.close_ledger(r)

    def test_stiff_small_capacity_limit(self):
        algebra=advance(S(2),P(0,10,30,50),100,.5,200,20)
        stiff=advance(S(2),P(1e-14,10,30,50),100,.5,200,20)
        self.assertAlmostEqual(stiff.state.temperature_C,algebra.state.temperature_C)
        self.assertAlmostEqual(stiff.total_ice_heat_J,algebra.total_ice_heat_J,places=9)
        self.close_ledger(stiff)

    def test_tiny_loss_retains_linear_heating_limit(self):
        r=advance(S(2),P(1000,1e-14,0,0),10,0,0,100)
        self.assertAlmostEqual(r.state.temperature_C,3)
        self.assertAlmostEqual(r.integrated_debris_temperature_K_s,25)
        self.close_ledger(r)

    def test_monotonic_subdivision_equivalence(self):
        p=P(3000,8,80,120)
        a=advance(S(3),p,1000,.2,400,10)
        state=S(3); ice=water=air=0
        for _ in range(100):
            b=advance(state,p,10,.2,400,10);state=b.state
            ice+=b.total_ice_heat_J;water+=b.water_to_debris_heat_J+b.water_to_bare_heat_J;air+=b.atmospheric_heat_J
        self.assertAlmostEqual(state.temperature_C,a.state.temperature_C,places=13)
        self.assertAlmostEqual(ice,a.total_ice_heat_J,places=8)
        self.assertAlmostEqual(water,a.water_to_debris_heat_J+a.water_to_bare_heat_J,places=8)

    def test_no_heat_or_storage_means_no_melt(self):
        for f in [0,1e-9,100,1e12]:
            r=advance(S(0),P(1000,1,100,100),50,0,f)
            self.assertEqual(r.total_ice_heat_J,0)
            self.close_ledger(r)

    def test_refreezing_and_cold_substrate_are_explicit(self):
        with self.assertRaises(UnsupportedThermalState):advance(S(0),P(10,2,0,0),10,0,0,-1)
        with self.assertRaises(UnsupportedThermalState):advance(S(1),P(10,2,0,0),10,0,0,0,ice_temperature_C=-1)
        with self.assertRaises(UnsupportedThermalState):advance(S(1),P(0,0,0,0),10,0,0,1)

    def test_bad_inputs_are_rejected(self):
        for value in [float('nan'),float('inf'),True,-1]:
            with self.assertRaises(ValueError):advance(S(0),P(1,1,1,1),10,0,value)
        with self.assertRaises(ValueError):advance(S(0),P(-1,1,1,1),10,0,1)
        with self.assertRaises(ValueError):advance(S(0),P(1,1,1,1),0,0,1)

    def test_independent_space_simpson_time_rk4(self):
        # Direct Simpson integration of the stream temperature profile, followed
        # by RK4 time integration of debris temperature and all flux inventories.
        # No analytic ODE/phi expression from the production kernel is used.
        C,K,Gb,Gd,F,Tin,A=3000.,9.,80.,120.,400.,.3,15.
        def rates(temp):
            star=Gd/(Gb+Gd)*temp
            qbare=qdebris=0.;N=64
            for j in range(N+1):
                x=j/N;tw=star+(Tin-star)*math.exp(-(Gb+Gd)*x/F)
                w=1 if j in [0,N] else (4 if j%2 else 2)
                qbare+=w*Gb*tw;qdebris+=w*Gd*(tw-temp)
            qb=qbare/(3*N);qd=qdebris/(3*N);qi=K*temp
            return [(A+qd-qi)/C,qb,qd,qi]
        y=[2.,0.,0.,0.];dt=.2
        for _ in range(1000):
            k1=rates(y[0]);k2=rates(y[0]+dt*k1[0]/2);k3=rates(y[0]+dt*k2[0]/2);k4=rates(y[0]+dt*k3[0])
            y=[v+dt*(a+2*b+2*c+d)/6 for v,a,b,c,d in zip(y,k1,k2,k3,k4)]
        r=advance(S(2),P(C,K,Gb,Gd),200,Tin,F,A)
        self.assertAlmostEqual(y[0],r.state.temperature_C,places=10)
        for independent,exact in [(y[1],r.bare_ice_heat_J),(y[2],r.water_to_debris_heat_J),(y[3],r.covered_ice_heat_J)]:
            self.assertLess(abs(independent-exact)/max(1,abs(exact)),5e-10)
        self.close_ledger(r)


if __name__=='__main__':unittest.main()
