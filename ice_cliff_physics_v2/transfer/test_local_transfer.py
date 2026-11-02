import math
import unittest
from local_transfer import ApplicabilityRefusal, ClosureConfig, config_from_dict, evaluate, finite_heat_chain, rectangular_flow


class GeometryTests(unittest.TestCase):
    def test_mass_and_reynolds_identity_for_different_aspect_ratios(self):
        for w,d,q in [(1.,.01,.003),(.2,.5,.004),(10.,.1,2.)]:
            r=rectangular_flow(q,w,d,1.8e-6)
            self.assertAlmostEqual(r['velocity_m_s']*r['area_m2'],q)
            self.assertAlmostEqual(r['Re'],r['velocity_m_s']*r['hydraulic_diameter_m']/1.8e-6)
            self.assertAlmostEqual(r['hydraulic_diameter_m']/4,r['area_m2']/r['wetted_perimeter_m'])

    def test_positive_flow_over_dry_geometry_is_not_silently_repaired(self):
        with self.assertRaises(ApplicabilityRefusal):rectangular_flow(.001,1.,0.)

    def test_zero_flow_is_not_natural_convection(self):
        for d in [0.,.1]:
            r=evaluate(0.,1.,d,{'method':'gnielinski_smooth'})
            self.assertEqual(r['h_W_m2_K'],0.)
            self.assertEqual(r['numerical_domain'],'zero_discharge_no_advective_exchange')
            self.assertEqual(r['standing_water_convection'],'not represented')

    def test_degenerate_and_nonfinite_inputs_refused(self):
        for q,w,d in [(True,1.,.1),(-1.,1.,.1),(1.,0.,.1),(1.,1.,float('nan')),(1.,1e308,1e308)]:
            with self.assertRaises(ApplicabilityRefusal):rectangular_flow(q,w,d)


class ClosureTests(unittest.TestCase):
    def test_legacy_fixed_coefficient_preserved(self):
        c=config_from_dict(None,legacy_fixed_h_W_m2_K=317.)
        self.assertEqual(evaluate(.01,1.,.1,c)['h_W_m2_K'],317.)
        self.assertEqual(evaluate(.01,1.,.1,{'method':'fixed'},legacy_fixed_h_W_m2_K=317.)['h_W_m2_K'],317.)

    def test_re_1e5_pr_one_reference_and_dimensions(self):
        # Q=.075, A=1, Dh=4/3 and nu=1e-6 => Re=1e5.
        # At Pr=1, Nu=.023*(1e5)**.8=230 exactly for either DB exponent.
        for method in ['db03','db04']:
            r=evaluate(.075,1.,1.,dict(method=method,nu_m2_s=1e-6,prandtl=1.,water_k_W_m_K=.5))
            self.assertAlmostEqual(r['Nu'],230.)
            self.assertAlmostEqual(r['h_W_m2_K'],86.25)

    def test_cooling_and_heating_are_distinct_choices(self):
        a=evaluate(.05,1.,.5,{'method':'db04'})
        b=evaluate(.05,1.,.5,{'method':'db03'})
        self.assertAlmostEqual(a['h_W_m2_K']/b['h_W_m2_K'],13.5**.1)

    def test_smooth_gnielinski_pr_one_simplification(self):
        r=evaluate(.075,1.,1.,dict(method='gnielinski_smooth',nu_m2_s=1e-6,prandtl=1.))
        self.assertTrue(.017 < r['smooth_tube_darcy_factor'] < .019)
        self.assertAlmostEqual(r['Nu'],r['smooth_tube_darcy_factor']/8*99000.)

    def test_low_reynolds_refused_without_hidden_laminar_switch(self):
        for method in ['db04','db03','gnielinski_smooth']:
            with self.assertRaises(ApplicabilityRefusal):evaluate(1e-6,1.,.1,dict(method=method))

    def test_high_re_needs_recorded_extrapolation(self):
        with self.assertRaises(ApplicabilityRefusal):evaluate(1.,1.,.5,dict(method='db03'))
        r=evaluate(1.,1.,.5,dict(method='db03',domain_policy='allow_high_re_extrapolation',extrapolation_reason='declared development sensitivity, not validation'))
        self.assertEqual(r['numerical_domain'],'explicit_high_Re_extrapolation')
        self.assertTrue(r['warnings'])

    def test_fitted_gain_needs_provenance_and_does_not_change_geometry(self):
        with self.assertRaises(ApplicabilityRefusal):ClosureConfig(gain=2.)
        a=evaluate(.05,1.,.5,dict(method='db03'))
        b=evaluate(.05,1.,.5,dict(method='db03',gain=2.,calibration_provenance='fixture-only training source'))
        self.assertEqual(a['Re'],b['Re'])
        self.assertEqual(b['h_W_m2_K'],2*a['h_W_m2_K'])

    def test_schema_unknown_fields_and_nonfinite_values_refused(self):
        for v in [dict(method='auto'),dict(schema_version=True),dict(fixed_h_W_m2_K=float('inf')),dict(target_h=3.),dict(domain_policy='allow_high_re_extrapolation')]:
            with self.assertRaises(ApplicabilityRefusal):config_from_dict(v)


class FiniteEnergyTests(unittest.TestCase):
    def test_serial_contact_equals_combined_contact_at_zero_surface_temperature(self):
        a=finite_heat_chain(.001,.5,[200.,300.,500.])
        b=finite_heat_chain(.001,.5,[1000.])
        self.assertAlmostEqual(a['delivered_heat_W'],b['delivered_heat_W'],places=12)
        self.assertAlmostEqual(a['outlet_heat_W'],b['outlet_heat_W'],places=12)

    def test_limiting_cases_and_global_energy_bound(self):
        for q,t,gs in [(0.,1.,[10.]),(.01,0.,[20.]),(.01,1.,[0.]),(1e-8,1.,[1e5]),(1e2,1.,[1e-8])]:
            r=finite_heat_chain(q,t,gs)
            self.assertLessEqual(r['delivered_heat_W'],r['inlet_heat_W']*(1+1e-14))
            self.assertGreaterEqual(r['outlet_heat_W'],0.)
            self.assertLessEqual(abs(r['residual_W']),1e-14*max(1.,r['inlet_heat_W']))

    def test_tiny_and_huge_ntu_retain_outgoing_heat(self):
        capacity=1000*4217*.001
        low=finite_heat_chain(.001,1.,[capacity*1e-12])
        high=finite_heat_chain(.001,1.,[capacity*50.])
        self.assertAlmostEqual(low['delivered_heat_W']/(capacity*1e-12),1.,places=10)
        self.assertGreater(high['outlet_heat_W'],0.)
        self.assertAlmostEqual(high['outlet_heat_W']/(capacity*math.exp(-50)),1.)

    def test_negative_or_nonfinite_conductance_refused(self):
        for gs in [[-1.],[float('nan')]]:
            with self.assertRaises(ApplicabilityRefusal):finite_heat_chain(.1,.1,gs)

    def test_finite_inputs_with_overflowing_derived_energy_are_refused(self):
        with self.assertRaises(ApplicabilityRefusal):finite_heat_chain(1e100,1e300,[1.])


if __name__=='__main__':unittest.main()
