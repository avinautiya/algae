import math
import unittest
import numpy as np
from dirichlet_slab import PreparedSlab,predict
from fv_reference import predict_fv


class SpectralTests(unittest.TestCase):
    def test_uniform_constant_is_exact(self):
        t=np.arange(51)*300.;z=np.array([.05,.1,.2,.3])
        r=predict(t,np.full(len(t),-5.),np.full(len(t),-5.),z[1:-1],3e-7,z,np.full(len(z),-5.))
        np.testing.assert_array_equal(r,np.full_like(r,-5.))

    def test_steady_linear_field_uses_real_endpoint_domain_not_ice_surface(self):
        t=np.arange(11)*300.;z=np.array([.05,.1,.15,.2]);initial=3+(z-.05)/.15*(-2)
        r=predict(t,np.full(len(t),3.),np.full(len(t),1.),z[1:-1],4e-7,z,initial)
        np.testing.assert_allclose(r,np.broadcast_to(initial[1:-1],r.shape),atol=2e-14,rtol=0)

    def test_independent_single_sine_decay(self):
        t=np.arange(21)*.1;z=np.linspace(0,1,513);target=np.array([.2,.5,.8]);initial=np.sin(math.pi*z)
        r=predict(t,np.zeros(len(t)),np.zeros(len(t)),target,.02,z,initial,modes=64,max_gap_s=.1)
        exact=np.exp(-.02*math.pi**2*t[:,None])*np.sin(math.pi*target[None,:])
        np.testing.assert_allclose(r[1:],exact[1:],atol=4e-6,rtol=0)

    def test_exact_ramp_particular_solution_uniform_and_irregular(self):
        alpha=.02;rate=.03;z=np.linspace(0,1,1025);initial=-rate*z*(1-z)/(2*alpha);target=np.array([.2,.5,.8])
        for t in [np.linspace(0,5,51),np.array([0.,.2,.5,1.,1.7,3.,5.])]:
            r=predict(t,rate*t,rate*t,target,alpha,z,initial,modes=64,max_gap_s=3.)
            exact=rate*t[:,None]-rate*target[None,:]*(1-target[None,:])/(2*alpha)
            np.testing.assert_allclose(r,exact,atol=3e-6,rtol=0)

    def test_independent_finite_volume_agrees_for_nontrivial_forcing(self):
        t=np.linspace(0,5,51);top=.8*np.sin(t/3);bottom=.2*np.sin(t/4)
        z=np.linspace(0,1,5);initial=np.array([0.,2.,-1.,1.,0.]);target=np.array([.2,.5,.8])
        r=predict(t,top,bottom,target,.02,z,initial,modes=128,max_gap_s=.1)
        f=predict_fv(t,top,bottom,target,.02,z,initial,cells=200,max_step_s=.0025)
        np.testing.assert_allclose(r[1:],f[1:],atol=6e-4,rtol=0)

    def test_depth_origin_and_temperature_offset_invariance(self):
        t=np.arange(12)*300.;z=np.array([.05,.1,.15,.2]);target=z[1:-1]
        top=np.sin(t/1000);bottom=.1*np.cos(t/2000);initial=np.array([top[0],-.5,.3,bottom[0]])
        a=predict(t,top,bottom,target,3e-7,z,initial)
        b=predict(t,top+100,bottom+100,target+10,3e-7,z+10,initial+100)
        np.testing.assert_allclose(a,b-100,atol=2e-11,rtol=0)

    def test_no_reset_or_future_interior_assimilation_contract(self):
        t=np.arange(10)*300.;z=np.array([0.,.1,.2]);top=np.arange(10,dtype=float);bottom=np.zeros(10)
        initial=np.array([0.,2.,0.]);prepared=PreparedSlab(t,top,bottom,[.1],z,initial)
        r=prepared.predict_with_diagnostics(1e-7)
        self.assertEqual(r.temperatures_C[0,0],2)
        self.assertIn('first input timestamp',r.diagnostics['initialization'])
        with self.assertRaises(ValueError):PreparedSlab(t,top,bottom,[.1],z,[top[-1],2.,0.])

    def test_gap_missing_and_boundary_target_refusals(self):
        z=[0.,.1,.2];t=[0.,300.,600.];top=[0.,1.,2.];bottom=[0.,0.,0.];initial=[0.,1.,0.]
        with self.assertRaises(ValueError):predict([0,300,901],top,bottom,[.1],1e-7,z,initial)
        with self.assertRaises(ValueError):predict(t,[0,np.nan,2],bottom,[.1],1e-7,z,initial)
        with self.assertRaises(ValueError):predict(t,top,bottom,[0.],1e-7,z,initial)
        with self.assertRaises(ValueError):predict([0,300,300],top,bottom,[.1],1e-7,z,initial)
        with self.assertRaises(ValueError):predict(t,top,bottom,[.1],-1e-7,z,initial)
        with self.assertRaises(ValueError):predict(t,top,bottom,[.1],1e-7,z,initial,modes=True)

    def test_preparation_does_not_alias_mutable_input_arrays(self):
        t=np.arange(10)*300.;top=np.linspace(0,1,10);bottom=np.zeros(10);z=np.array([0.,.1,.2]);initial=np.array([0.,1.,0.])
        p=PreparedSlab(t,top,bottom,[.1],z,initial);a=p.predict(1e-7)
        top[:]=99;initial[:]=99;t[:]=0
        np.testing.assert_array_equal(p.predict(1e-7),a)

if __name__=='__main__':unittest.main()
