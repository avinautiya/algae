import math,unittest
from thermal_identifiability import harmonic_temperature,harmonic_bottom_flux,harmonic_weighted_moment
class IdentifiabilityTests(unittest.TestCase):
    def test_boundaries(self):
        for t in [0,3600,12345]:
            self.assertAlmostEqual(harmonic_temperature(0,t,length_m=.2,alpha_m2_s=1e-6),5+2*math.cos(2*math.pi*t/86400),places=12)
            self.assertAlmostEqual(harmonic_temperature(.2,t,length_m=.2,alpha_m2_s=1e-6),0,places=12)
    def test_k_and_capacity_scaling_not_identified_by_temperature(self):
        k,C=1.,1e6;scale=7.;alpha=k/C
        a=harmonic_temperature(.1,43200,length_m=.2,alpha_m2_s=alpha)
        b=harmonic_temperature(.1,43200,length_m=.2,alpha_m2_s=scale*k/(scale*C))
        self.assertEqual(a,b)
        q1=harmonic_bottom_flux(43200,length_m=.2,alpha_m2_s=alpha,conductivity_W_m_K=k)
        q2=harmonic_bottom_flux(43200,length_m=.2,alpha_m2_s=alpha,conductivity_W_m_K=scale*k)
        self.assertAlmostEqual(q2,scale*q1,places=12)
    def test_local_moment_identity(self):
        L=.2;alpha=1e-6;k=1.3;Cv=k/alpha
        for t in [1000.,10000.,36000.]:
            eps=.1
            deriv=(harmonic_weighted_moment(t+eps,length_m=L,alpha_m2_s=alpha)-harmonic_weighted_moment(t-eps,length_m=L,alpha_m2_s=alpha))/(2*eps)
            qs=k*(5+2*math.cos(2*math.pi*t/86400))/L
            qb=harmonic_bottom_flux(t,length_m=L,alpha_m2_s=alpha,conductivity_W_m_K=k)
            self.assertAlmostEqual(qs-Cv*deriv,qb,places=8)
    def test_period_integrated_moment_cancels(self):
        a=harmonic_weighted_moment(0,length_m=.2,alpha_m2_s=1e-6)
        b=harmonic_weighted_moment(86400,length_m=.2,alpha_m2_s=1e-6)
        self.assertAlmostEqual(a,b,places=14)
    def test_invalid(self):
        with self.assertRaises(ValueError):harmonic_temperature(-.1,0,length_m=.2,alpha_m2_s=1e-6)
        with self.assertRaises(ValueError):harmonic_temperature(.1,0,length_m=.2,alpha_m2_s=0)
if __name__=='__main__':unittest.main()
