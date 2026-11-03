import math,random,unittest
from ice_reference import impedance,identify,enclosure
def context():return {k:True for k in ['entire_transformed_history_nonmelting','initial_ice_temperature_anomaly_zero','bottom_temperature_anomaly_zero','common_clock_and_transform_convention','debris_boundary_ratio_at_known_ice_interface','ice_properties_independent_of_melt_target']}
KW=dict(ice_conductivity_W_mK=2.1,ice_capacity_J_m3K=917*2100.,ice_depth_m=2.)
class ReferenceTests(unittest.TestCase):
    def test_reconstruct_known_parameters(self):
        s1,s2=1/86400.,1/(4*86400.)
        for k,C,R in [(1.2,1.8e6,0.),(2.,2.2e6,.03),(.3,1e6,.2)]:
            a1=k*(R+impedance(s1,**KW));a2=k*(R+impedance(s2,**KW))
            r=identify(a1,a2,s1,s2,diffusivity_m2_s=k/C,context=context(),**KW)
            self.assertAlmostEqual(r['conductivity_W_mK'],k,places=12)
            self.assertAlmostEqual(r['volumetric_capacity_J_m3K']/C,1,places=12)
            self.assertAlmostEqual(r['contact_resistance_m2K_W'],R,places=12)
    def test_finite_depth_limits(self):
        self.assertAlmostEqual(impedance(1e-16,**KW),2/2.1,places=8)
        r=impedance(1.,**KW);half=impedance(1.,**dict(KW,ice_depth_m=math.inf))
        self.assertAlmostEqual(r,half,places=15)
    def test_unresolved_contrast_not_repaired(self):
        z=[impedance(s,**KW)for s in[1/86400.,1/(4*86400.)]]
        r=enclosure([.05,.5],[.05,.5],1/86400.,1/(4*86400.),diffusivity_interval_m2_s=[1e-7,1e-6],context=context(),**KW)
        self.assertEqual(r['status'],'positive_conductivity_not_resolved_by_bounds')
    def test_conditional_box_contains_true_samples(self):
        random.seed(20261103);s1=1/(4*86400.);s2=1/86400.;k=1.2;C=1.8e6;R=.03
        a1=k*(R+impedance(s1,**KW));a2=k*(R+impedance(s2,**KW));eps=.002
        e=enclosure([a1-eps,a1+eps],[a2-eps,a2+eps],s1,s2,diffusivity_interval_m2_s=[.99*k/C,1.01*k/C],context=context(),**KW)
        self.assertEqual(e['status'],'conditionally_enclosed')
        for _ in range(100):
            r=identify(a1+random.uniform(-eps,eps),a2+random.uniform(-eps,eps),s1,s2,diffusivity_m2_s=random.uniform(.99*k/C,1.01*k/C),context=context(),**KW)
            for field,box in [('conductivity_W_mK','conductivity_W_mK'),('volumetric_capacity_J_m3K','capacity_J_m3K'),('contact_resistance_m2K_W','contact_resistance_m2K_W')]:
                self.assertLessEqual(e[box][0],r[field]);self.assertLessEqual(r[field],e[box][1])
    def test_warm_scaling_does_not_resolve_absolute_scale(self):
        k,C,R=1.2,1.8e6,.03
        for scale in [.5,2,5]:
            self.assertAlmostEqual(k*scale/(C*scale),k/C,places=18)
            self.assertAlmostEqual(k*scale*(R/scale),k*R)
    def test_missing_history_and_equal_rates_refused(self):
        c=context();c['entire_transformed_history_nonmelting']=False
        with self.assertRaises(ValueError):identify(.2,.3,1e-5,2e-5,diffusivity_m2_s=1e-6,context=c,**KW)
        with self.assertRaises(ValueError):identify(.2,.3,1e-5,1e-5,diffusivity_m2_s=1e-6,context=context(),**KW)
if __name__=='__main__':unittest.main()
