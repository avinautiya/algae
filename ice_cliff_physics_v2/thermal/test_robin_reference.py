import math
import unittest
from robin_reference import RobinSlab,eigenvalues,water_exchange_conductance
from resolved_slab import SlabState,SlabParameters,advance


class RobinTests(unittest.TestCase):
    def test_eigenvalues_obey_intervals_and_robin_equation(self):
        for B in [0.,.1,1.,10.,1e6]:
            roots=eigenvalues(B,32)
            for n,x in enumerate(roots):
                self.assertGreaterEqual(x,(n+.5)*math.pi)
                self.assertLessEqual(x,(n+1)*math.pi)
                self.assertLess(abs(x*math.cos(x)+B*math.sin(x))/(1+x+B),2e-13)

    def test_initial_and_steady_limits(self):
        a=RobinSlab(1.,1.,1.,.5,terms=128);b=RobinSlab(1.,1.,1.,.5,terms=512)
        self.assertEqual(a.response(0)['stored_energy_J'],0)
        self.assertEqual(a.temperature(.5,0),0)
        # Finite-mode initial surface residual is exposed, not hidden by the
        # exact initial-state special case. It converges only algebraically.
        self.assertLess(abs(b.initial_series_residual(1)),abs(a.initial_series_residual(1))/3.9)
        for x in [0.,.2,.5,.8,1.]:self.assertAlmostEqual(b.temperature(x,10),b.Ts*x,places=10)
        self.assertAlmostEqual(b.response(10)['ice_power_W'],b.K*b.Ts,places=10)

    def test_boundary_and_independent_finite_difference_pde(self):
        s=RobinSlab.from_contact(90000,20,400,200,.5,3,512)
        for fo in [.03,.2,1.]:
            t=fo*s.tau
            self.assertEqual(s.temperature(0,t),0)
            boundary=s.K*s.spatial_derivative(1,t)-(s.P-s.H*s.temperature(1,t))
            self.assertLess(abs(boundary),1e-10)
            for x in [.2,.5,.8]:
                dt=1e-5*s.tau;dx=1e-4
                derivative=(s.temperature(x,t+dt)-s.temperature(x,t-dt))/(2*dt)
                curvature=(s.temperature(x+dx,t)-2*s.temperature(x,t)+s.temperature(x-dx,t))/dx**2
                self.assertLess(abs(derivative-curvature/s.tau)/(s.Ts/s.tau),5e-6)

    def test_energy_and_256_512_mode_stability(self):
        for qa in [0.,20.,-10.]:
            a=RobinSlab.from_contact(90000,20,400,200,.5,qa,256)
            b=RobinSlab.from_contact(90000,20,400,200,.5,qa,512)
            for fo in [.01,.1,1.,10.]:
                x=a.response(fo*a.tau);y=b.response(fo*b.tau)
                scale=max(1.,abs(y['boundary_heat_J']),abs(y['stored_energy_J']),abs(y['ice_heat_J']))
                self.assertLess(abs(y['energy_residual_J'])/scale,1e-12)
                self.assertLess(abs(x['ice_heat_J']-y['ice_heat_J'])/scale,1e-12)
                self.assertLess(abs(x['stored_energy_J']-y['stored_energy_J'])/scale,1e-12)

    def test_fv_time_converges_to_analytic_reference(self):
        p=SlabParameters(1,1,1,1);s=RobinSlab.from_contact(1,1,10,2,.5,0,512);errors=[]
        for steps in [20,40,80]:
            state=SlabState(tuple([0.]*64));heat=0.
            for _ in range(steps):
                r=advance(state,p,.2/steps,.5,10,debris_water_conductance_W_K=2)
                state=r['state'];heat+=r['covered_ice_heat_J']
            errors.append(abs(heat-s.response(.2)['ice_heat_J']))
        self.assertGreater(errors[0]/errors[1],1.8)
        self.assertGreater(errors[1]/errors[2],1.8)

    def test_fv_mesh_converges_to_independent_cell_averages(self):
        p=SlabParameters(1,1,1,1);s=RobinSlab.from_contact(1,1,10,2,.5,0,512);errors=[]
        for n in [8,16,32]:
            state=SlabState(tuple([0.]*n))
            for _ in range(4000):state=advance(state,p,.2/4000,.5,10,debris_water_conductance_W_K=2)['state']
            exact=[s.cell_average(i/n,(i+1)/n,.2) for i in range(n)]
            errors.append(math.sqrt(math.fsum((a-b)**2 for a,b in zip(state.temperatures_C,exact))/n))
        self.assertGreater(errors[0]/errors[1],3)
        self.assertGreater(errors[1]/errors[2],2.5)

    def test_generic_boundary_handles_parallel_bare_contact_without_false_water_split(self):
        F,Gb,Gd,Tin,qa=200.,50.,80.,.5,3.
        G=Gb+Gd;phi=-math.expm1(-G/F)/(G/F);alpha=Gd/G
        a=Gd*phi;b=Gd*(Gb/G+phi*alpha)
        s=RobinSlab.from_linear_boundary(90000,20,qa+a*Tin,b,512)
        p=SlabParameters(.05,1,1.8e6,1);state=SlabState(tuple([0.]*128));heat=0.
        for _ in range(800):
            r=advance(state,p,s.tau*.2/800,Tin,F,Gb,Gd,qa);state=r['state'];heat+=r['covered_ice_heat_J']
        exact=s.response(.2*s.tau)
        self.assertLess(abs(heat-exact['ice_heat_J'])/exact['ice_heat_J'],.01)
        self.assertIsNone(exact['water_heat_J'])
        self.assertLess(abs(exact['energy_residual_J']),1e-8)

    def test_dry_neumann_limit(self):
        s=RobinSlab.from_contact(1,1,0,200,0,2,512)
        self.assertEqual(s.H,0)
        self.assertEqual(s.roots[0],math.pi/2)
        out=s.response(10)
        self.assertEqual(out['water_heat_J'],0)
        self.assertAlmostEqual(out['ice_power_W'],2,places=9)
        self.assertAlmostEqual(out['boundary_heat_J'],20)

    def test_no_forcing_stays_zero(self):
        s=RobinSlab(100,20,30,0,0,256)
        for t in [0.,1.,1000.]:
            r=s.response(t)
            for k in ['ice_power_W','ice_heat_J','stored_energy_J','boundary_heat_J','energy_residual_J']:self.assertEqual(r[k],0)

    def test_bad_inputs_refuse(self):
        for args in [(0,1,1,.5),(1,-1,1,.5),(1,1,-1,.5),(1,1,1,-.5)]:
            with self.assertRaises(ValueError):RobinSlab(*args)
        with self.assertRaises(ValueError):RobinSlab(1,1,1,.5,-1)
        with self.assertRaises(ValueError):eigenvalues(1,0)
        with self.assertRaises(ValueError):RobinSlab(1,1,1,.5).response(-1)
        with self.assertRaises(ValueError):RobinSlab(1,1,1,.5).temperature(2,1)
        self.assertEqual(water_exchange_conductance(0,100),0)

if __name__=='__main__':unittest.main()
