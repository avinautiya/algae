import math
import unittest
from resolved_slab import SlabState,SlabParameters,ArealSlabParameters,advance


class SlabTests(unittest.TestCase):
    def test_steady_robin_solution_is_exact_on_each_mesh(self):
        p=SlabParameters(.05,1.,1.8e6,1.)
        F,G,Tw=400.,200.,.5
        H=F*(-math.expm1(-G/F));surface=H*Tw/(H+20)
        for n in [1,2,8,32]:
            state=SlabState(tuple(surface*(i+.5)/n for i in range(n)))
            r=advance(state,p,1000,Tw,F,debris_water_conductance_W_K=G)
            self.assertLess(max(abs(a-b) for a,b in zip(r['state'].temperatures_C,state.temperatures_C)),1e-13)
            self.assertAlmostEqual(r['surface_temperature_C'],surface,places=12)
            self.assertAlmostEqual(r['covered_ice_heat_J']/1000,20*surface,places=11)
            self.assertLess(abs(r['global_energy_residual_J']),1e-8)

    def test_closed_discrete_energy_with_parallel_contacts(self):
        p=SlabParameters(.03,.8,2e6,2.)
        state=SlabState(tuple([4.]*17))
        for _ in range(30):
            r=advance(state,p,13,.1,100,80,200,40);state=r['state']
            scale=max(1.,abs(r['debris_storage_change_J']),r['total_ice_heat_J'])
            self.assertLess(abs(r['global_energy_residual_J'])/scale,2e-12)
            self.assertLess(abs(r['water_exchange_residual_J'])/scale,2e-12)
            self.assertGreaterEqual(r['water_outlet_heat_J'],0)

    def test_dry_limit_has_no_water_heat(self):
        r=advance(SlabState((2.,2.,2.)),SlabParameters(.05,1,1e6,1),10,0,0,100,100,20)
        self.assertAlmostEqual(r['water_to_debris_heat_J'],0)
        self.assertEqual(r['bare_ice_heat_J'],0)
        self.assertEqual(r['water_inlet_heat_J'],0)
        self.assertLess(abs(r['global_energy_residual_J']),1e-8)

    def test_analytic_diffusion_limit_mesh_converges(self):
        # Unit slab initially zero, bottom T=0/top T=1. The cell-average exact
        # sine-series solution supplies an independent PDE reference.
        errors=[];t=.1
        for n in [8,16,32]:
            p=SlabParameters(1,1,1,1);state=SlabState(tuple([0.]*n))
            for _ in range(2000):
                state=advance(state,p,t/2000,1,1e12,debris_water_conductance_W_K=1e8)['state']
            exact=[]
            for i in range(n):
                lo=i/n;hi=(i+1)/n;mean=(lo+hi)/2
                for k in range(1,30):
                    mean+=2*(-1)**k/(k*math.pi)*math.exp(-(k*math.pi)**2*t)*(math.cos(k*math.pi*lo)-math.cos(k*math.pi*hi))/(k*math.pi*(hi-lo))
                exact.append(mean)
            errors.append(math.sqrt(sum((a-b)**2 for a,b in zip(state.temperatures_C,exact))/n))
        self.assertGreater(errors[0]/errors[1],3)
        self.assertGreater(errors[1]/errors[2],2.5)

    def test_areal_coordinates_match_physical_slab_without_new_parameter(self):
        physical=SlabParameters(.05,1.,1.8e6,2.)
        areal=ArealSlabParameters(90000.,20.,2.)
        state=SlabState(tuple([1.]*12))
        a=advance(state,physical,100,.3,200,60,80,20)
        b=advance(state,areal,100,.3,200,60,80,20)
        self.assertEqual(a['state'],b['state'])
        self.assertEqual(a['covered_ice_heat_J'],b['covered_ice_heat_J'])
        self.assertAlmostEqual(physical.diffusion_time_s,areal.diffusion_time_s,places=10)
        self.assertIsNone(b['cell_width_m'])

    def test_hot_dry_storage_release_matches_analytic_series(self):
        n=64;state=SlabState(tuple([4.]*n));p=SlabParameters(1,1,1,1);heat=0.
        for _ in range(2000):
            r=advance(state,p,.1/2000,0,0);state=r['state'];heat+=r['covered_ice_heat_J']
        remaining=sum(8/(math.pi**2*(2*j+1)**2)*math.exp(-((2*j+1)*math.pi)**2*.1/4) for j in range(100))
        expected=4*(1-remaining)
        self.assertLess(abs(heat-expected)/expected,5e-4)

    def test_time_refinement_is_first_order_for_diffusion_eigenmode(self):
        n=64;p=SlabParameters(1,1,1,1);initial=tuple(math.sin(math.pi*(i+.5)/(2*n)) for i in range(n));errors=[]
        exact=[v*math.exp(-math.pi**2*.2/4) for v in initial]
        for steps in [10,20,40]:
            state=SlabState(initial)
            for _ in range(steps):state=advance(state,p,.2/steps,0,0)['state']
            errors.append(max(abs(a-b) for a,b in zip(state.temperatures_C,exact)))
        self.assertGreater(errors[0]/errors[1],1.8)
        self.assertGreater(errors[1]/errors[2],1.8)

    def test_invalid_states_refuse(self):
        with self.assertRaises(ValueError):advance(SlabState((-1.,)),SlabParameters(.05,1,1e6,1),10,0,1)
        with self.assertRaises(ValueError):advance(SlabState((1.,)),SlabParameters(.05,1,1e6,1),0,0,1)
        with self.assertRaises(ValueError):advance(SlabState((0.,)),SlabParameters(.05,1,1e6,1),10,0,0,atmospheric_power_W=-1)


if __name__=='__main__':unittest.main()
