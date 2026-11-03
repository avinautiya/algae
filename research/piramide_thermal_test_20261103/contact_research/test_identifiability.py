import unittest
from contact_column import simulate
from identifiability import invariants,steady_temperate_sensor_bound


class ContactIdentifiabilityTests(unittest.TestCase):
    def test_exact_warm_temperature_scaling_invariance_but_melt_changes(self):
        k=.7;C=1e6;R=.2;scale=3.2
        def run(k,C,R):
            return simulate(debris_conductivity_W_mK=k,debris_heat_capacity_J_m3K=C,debris_depth_m=.1,
                initial_debris_temperature_C=[0.]*12,initial_ice_temperature_C=[0.]*24,
                surface_temperature_C=[4.,10.,2.],interval_seconds=900.,dt_seconds=30.,
                ice_depth_m=.2,contact_resistance_m2K_W=R,record_profiles=True)
        a=run(k,C,R);b=run(k*scale,C*scale,R/scale)
        for x,y in zip(a['history'],b['history']):
            self.assertEqual(x['interface_temperature_C'],0.)
            self.assertEqual(y['interface_temperature_C'],0.)
            self.assertLess(max(abs(u-v) for u,v in zip(x['debris_temperature_C'],y['debris_temperature_C'])),3e-13)
        self.assertAlmostEqual(b['melt_mmwe']/a['melt_mmwe'],scale,places=12)
        aa=invariants(k,C,R,.1);bb=invariants(k*scale,C*scale,R/scale,.1)
        self.assertAlmostEqual(aa['diffusivity_m2_s'],bb['diffusivity_m2_s'])
        self.assertAlmostEqual(aa['robin_length_m'],bb['robin_length_m'])

    def test_steady_geometry_contact_degeneracy(self):
        # Both geometries lie within a nominal .50m +/- .02m thickness range.
        a=invariants(1.,2e6,.1,.5);b=invariants(1.,2e6,.08,.52)
        self.assertAlmostEqual(a['steady_effective_depth_m'],b['steady_effective_depth_m'])
        for z in (.125,.25,.375,.48):
            self.assertAlmostEqual(10*(1-z/a['steady_effective_depth_m']),10*(1-z/b['steady_effective_depth_m']))

    def test_conditional_sensor_bound_and_no_unsupported_phase_claim(self):
        r=steady_temperate_sensor_bound(sensor_temperature_C=.374,sensor_absolute_bound_C=.4,
            mean_heat_flux_W_m2=23.5,temperate_steady_assumption_explicit=True)
        self.assertAlmostEqual(r['necessary_R_upper_bound_m2K_W'],.774/23.5)
        self.assertLess(r['necessary_R_upper_bound_m2K_W'],.371)
        with self.assertRaises(ValueError):steady_temperate_sensor_bound(sensor_temperature_C=.374,sensor_absolute_bound_C=.4,
            mean_heat_flux_W_m2=23.5,temperate_steady_assumption_explicit=False)


if __name__=='__main__':unittest.main()
