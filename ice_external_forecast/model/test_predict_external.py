import copy
import datetime as dt
import math
import unittest
from predict_external import (Section,predict_interval,predict_packet,ApplicabilityRefusal,
    PHYSICAL_MODELS,REFERENCE,RHO_W,CP_W,RHO_I,LATENT,CP_I,H_FIXED,RATE)

def reach(length=25):
    return Section('synthetic reach',2.,.5,length,False,'synthetic fixture','synthetic_fixture')

def packet():
    return dict(purpose='synthetic_verification',clock_convention='UTC',
        alignment=dict(path_alignment_verified=True,geometry_independent_of_bed_target=True,source='synthetic fixture'),
        records=[dict(start='2000-01-01T00:00:00Z',end='2000-01-01T01:00:00Z',Q_m3_s=.5,
            inlet_temperature_C=1.,ice_temperature_C=0.,forcing_source='synthetic fixture',forcing_basis='synthetic_fixture',
            sections=[dict(name='reach',width_m=2.,depth_m=.5,length_m=25.,source='synthetic fixture',geometry_basis='synthetic_fixture')])])

class ExternalAdapterChecks(unittest.TestCase):
    def test_uniform_reach_exact_exponential(self):
        s=reach();r=predict_interval([s],.5,1.,3600.,0.,'constant_h_fitted')
        expected=math.exp(-H_FIXED*3.*25./(RHO_W*CP_W*.5))
        self.assertAlmostEqual(r['sections'][0]['outlet_temperature_C'],expected)
        self.assertTrue(math.isclose(r['removed_water_heat_J'],RHO_W*CP_W*.5*(1-expected)*3600.,rel_tol=1e-14))

    def test_splitting_identical_geometry_is_invariant(self):
        for model in PHYSICAL_MODELS:
            one=predict_interval([reach(50)],.5,.7,100.,0.,model)
            two=predict_interval([reach(20),reach(30)],.5,.7,100.,0.,model)
            self.assertAlmostEqual(one['outlet_heat_J'],two['outlet_heat_J'],places=6)
            self.assertAlmostEqual(one['removed_ice_volume_m3'],two['removed_ice_volume_m3'],places=12)

    def test_cold_ice_latent_sensible_split_and_footprint(self):
        warm=predict_interval([reach()],.5,1.,100.,0.,'constant_h_fitted')
        cold=predict_interval([reach()],.5,1.,100.,-5.,'constant_h_fitted')
        self.assertAlmostEqual(cold['removed_ice_volume_m3']/warm['removed_ice_volume_m3'],LATENT/(LATENT+5*CP_I))
        row=cold['sections'][0]
        self.assertAlmostEqual(row['bed_heat_J']/row['delivered_heat_J'],2/3)
        self.assertAlmostEqual(row['sidewall_ice_volume_m3']+row['bed_ice_volume_m3'],row['removed_ice_volume_m3'])
        self.assertAlmostEqual(row['mean_bed_lowering_m']*2*25,row['bed_ice_volume_m3'])
        self.assertLess(abs(cold['interface_energy_residual_J'])/cold['removed_water_heat_J'],1e-14)

    def test_zero_inlet_heat_all_physical_models(self):
        for model in PHYSICAL_MODELS:
            r=predict_interval([reach()],.5,0.,3600.,0.,model)
            self.assertEqual(r['removed_ice_volume_m3'],0.)
            self.assertEqual(r['sections'][0]['mean_bed_lowering_m'],0.)

    def test_large_ntu_cannot_use_more_than_inlet_heat(self):
        for model in PHYSICAL_MODELS:
            r=predict_interval([reach(1e6)],.5,1.,3600.,0.,model)
            self.assertLessEqual(r['removed_water_heat_J'],r['inlet_heat_J']*(1+1e-14))
            self.assertGreaterEqual(r['outlet_heat_J'],0.)
            self.assertLess(abs(r['energy_residual_J'])/r['inlet_heat_J'],1e-14)

    def test_tiny_ntu_uses_stable_expm1(self):
        r=predict_interval([reach(1e-12)],.5,1.,1.,0.,'constant_h_fitted')
        self.assertGreater(r['removed_water_heat_J'],0.)
        self.assertAlmostEqual(r['removed_water_heat_J'],H_FIXED*3e-12,places=16)

    def test_critical_point_is_local_operator_not_finite_volume(self):
        p=Section.from_input(dict(name='control',kind='point_control',length_m=0.,width_m=2.,
            depth_mode='critical_control',critical_control_verified=True,critical_control_source='synthetic control',
            source='synthetic fixture',geometry_basis='synthetic_fixture'),.5)
        r=predict_interval([reach(),p],.5,1.,100.,0.,'constant_h_fitted')
        point=r['sections'][1]
        self.assertAlmostEqual(point['Fr'],1.)
        self.assertEqual(point['removed_ice_volume_m3'],0.)
        self.assertEqual(point['delivered_heat_J'],0.)
        self.assertGreater(point['point_bed_lowering_m'],0.)
        self.assertEqual(point['inlet_temperature_C'],r['sections'][0]['outlet_temperature_C'])

    def test_missing_field_route_and_critical_evidence_refused(self):
        p=packet();p['alignment']['path_alignment_verified']=False
        with self.assertRaises(ApplicabilityRefusal):predict_packet(p)
        p=packet();del p['records'][0]['sections'][0]['depth_m']
        with self.assertRaises(ApplicabilityRefusal):predict_packet(p)
        with self.assertRaises(ApplicabilityRefusal):Section.from_input(dict(name='bad',kind='uniform_reach',length_m=20,width_m=2,depth_mode='critical_control',source='synthetic',geometry_basis='synthetic_fixture',critical_control_verified=True),.5)

    def test_time_gap_and_clock_invention_refused(self):
        p=packet();second=copy.deepcopy(p['records'][0]);second.update(start='2000-01-01T02:00:00Z',end='2000-01-01T03:00:00Z');p['records'].append(second)
        with self.assertRaises(ApplicabilityRefusal):predict_packet(p)
        p=packet();p['records'][0]['start']='2000-01-01T00:00:00'
        with self.assertRaises(ApplicabilityRefusal):predict_packet(p)

    def test_uncapped_reference_identity_and_energy_warning_preserved(self):
        p=packet();p['records'][0]['inlet_temperature_C']=0.
        out=predict_packet(p);models=out['intervals'][0]['models'];self.assertEqual(len(models),7)
        ref=models[REFERENCE]
        self.assertAlmostEqual(ref['sections'][0]['normal_or_vertical_bed_lowering_m'],RATE*3600.)
        self.assertTrue(ref['exceeds_inlet_heat_supply']);self.assertFalse(ref['closed_physical_component_ledger'])
        self.assertEqual(models['constant_rate_fitted_heat_limited']['removed_ice_volume_m3'],0.)

    def test_gravitational_scale_does_not_add_unmeasured_heat(self):
        baseline=predict_interval([reach()],.5,.7,100.,0.,'constant_h_fitted')
        s=Section('potential example',2.,.5,25.,False,'synthetic fixture','synthetic_fixture',3.)
        diagnostic=predict_interval([s],.5,.7,100.,0.,'constant_h_fitted')
        self.assertEqual(baseline['removed_ice_volume_m3'],diagnostic['removed_ice_volume_m3'])
        self.assertGreater(diagnostic['sections'][0]['potential_energy_power_scale_W'],0.)

    def test_synthetic_forcing_and_invalid_section_cannot_be_promoted(self):
        p=packet();p['purpose']='conditional_external_hindcast'
        with self.assertRaises(ApplicabilityRefusal):predict_packet(p)
        with self.assertRaises(ApplicabilityRefusal):reach(-1)

    def test_finite_inputs_with_underflowing_derived_geometry_refused(self):
        tiny=Section('tiny',1e-300,1e-300,1.,False,'synthetic','synthetic_fixture')
        with self.assertRaises(ApplicabilityRefusal):tiny.geometry(.5)

if __name__=='__main__':unittest.main()
