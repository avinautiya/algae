import unittest
import numpy as np
from dirichlet_slab import PreparedSlab

class BottomGradientTests(unittest.TestCase):
    def test_steady_signed(self):
        t=np.arange(6.)*100; z=np.linspace(0,.2,51)
        for top,bot in [(3.,1.),(-1.,2.)]:
            p=PreparedSlab(t,np.full(6,top),np.full(6,bot),[.1],z,np.linspace(top,bot,51))
            np.testing.assert_allclose(p.bottom_downward_gradient(1e-6).gradient_K_m,(top-bot)/.2,atol=1e-12)
            np.testing.assert_allclose(p.bottom_interval_gradient(1e-6).integral_K_s_m,(top-bot)/.2*100,atol=1e-10)
    def test_sine_decay_and_moment(self):
        L=.2; alpha=1e-6;t=np.arange(101.)*100; z=np.linspace(0,L,20001)
        p=PreparedSlab(t,np.zeros(len(t)),np.zeros(len(t)),[.1],z,np.sin(np.pi*z/L),modes=128)
        rate=alpha*(np.pi/L)**2
        expected=np.pi/L*np.exp(-rate*t)
        np.testing.assert_allclose(p.bottom_downward_gradient(alpha).gradient_K_m,expected,rtol=1e-7)
        expectedI=np.pi/L/rate*(np.exp(-rate*t[:-1])-np.exp(-rate*t[1:]))
        np.testing.assert_allclose(p.bottom_interval_gradient(alpha).integral_K_s_m,expectedI,rtol=1e-7)
    def test_ramp_convergence(self):
        L=.2; alpha=1e-6;r=1e-4;t=np.arange(101.)*100;z=np.linspace(0,L,20001)
        initial=-r*z*(L-z)/(2*alpha);expected=-r*L/(2*alpha)
        errors=[]
        for n in [32,64,128]:
            p=PreparedSlab(t,r*t,r*t,[.1],z,initial,modes=n)
            errors.append(abs(p.bottom_downward_gradient(alpha).gradient_K_m[-1]-expected))
            np.testing.assert_allclose(p.bottom_interval_gradient(alpha).interval_mean_K_m,expected,rtol=2e-5)
        self.assertGreater(errors[0]/errors[1],1.9);self.assertGreater(errors[1]/errors[2],1.9)
    def test_assumed_provenance(self):
        p=PreparedSlab.with_assumed_ice_bottom([0,300],[2,2],[.05],.3,[0,.1],[2,1],'synthetic independent depth')
        d=p.predict_with_diagnostics(1e-6).diagnostics
        self.assertIsNone(d['measured_domain_depths_m']);self.assertEqual(d['domain_provenance']['bottom_boundary'],'assumed_ice_Tm0')
        with self.assertRaises(ValueError):PreparedSlab.with_assumed_ice_bottom([0,300],[2,2],[.05],.1,[0,.1],[2,1],'depth')

if __name__=='__main__':unittest.main()
