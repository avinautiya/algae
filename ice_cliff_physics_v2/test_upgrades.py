import copy
import json
import unittest
from pathlib import Path
from model import forecast


def warm_case():
    c=json.loads((Path(__file__).parent/'example.json').read_text())
    c['physics'].update(initial_surface_cold_J_m2=0.,new_surface_cold_J_m2=0.,ice_subcooling_K=0.,face_cover_fraction=.02)
    c['periods']=[dict(duration_s=600.,inflow_m3_s=.07,water_superheat_K=.5,
        face_atmospheric_net_W_m2=100.,floor_atmospheric_net_W_m2=100.,bypass_atmospheric_net_W_m2=100.)]
    c['step_s']=60.;c['removal_time_s']=None
    c['upgrades']=dict(wet_debris='resolved_temperate',face_cover_location='lower')
    return c


class UpgradeTests(unittest.TestCase):
    def test_lower_cover_wets_without_flooding_tall_cliff(self):
        r=forecast(warm_case())
        self.assertTrue(all(x['wet_covered_face_area_m2']>0 for x in r['series']))
        self.assertTrue(all(x['contact_depth_m']<.5 for x in r['series']))
        self.assertLess(abs(r['checks']['global_energy_residual_J']),1e-5)
        self.assertLess(abs(r['checks']['geometry_phase_volume_residual_m3']),1e-12)

    def test_partial_wetting_does_not_mix_entire_debris_temperature(self):
        c=warm_case();c['physics']['face_cover_fraction']=.2
        r=forecast(c);p=r['series'][-1]['wet_debris_thermal']['partition']
        self.assertGreater(p['wet_area_m2'],0);self.assertGreater(p['dry_area_m2'],0)
        self.assertNotEqual(p['wet_profile_C'],p['dry_profile_C'])

    def test_cold_wet_ice_is_rejected_before_running(self):
        c=warm_case();c['physics']['ice_subcooling_K']=1.
        with self.assertRaisesRegex(ValueError,'temperate'):forecast(c)

    def test_zero_inflow_solar_melt_has_explicit_downstream_mass(self):
        c=warm_case();c['physics']['face_cover_fraction']=0.;c['periods'][0]['inflow_m3_s']=0.
        c['upgrades']['meltwater_policy']=dict(mode='downstream_free_drainage',max_generated_discharge_m3_s=.001,
            provenance='SYNTHETIC unponded downstream drainage capacity; not a solved flow field')
        r=forecast(c)
        self.assertGreater(r['totals']['phase_water_m3'],0.)
        self.assertEqual(r['totals']['water_in_m3'],0.)
        self.assertAlmostEqual(r['totals']['water_out_m3'],r['totals']['phase_water_m3'])
        self.assertTrue(all(x['meltwater_to_inflow_ratio'] is None for x in r['series']))

    def test_zero_inflow_requires_capacity_not_an_infinite_ratio_hack(self):
        c=warm_case();c['physics']['face_cover_fraction']=0.;c['periods'][0]['inflow_m3_s']=0.
        with self.assertRaisesRegex(ValueError,'free-drainage'):forecast(c)
        c['upgrades']['meltwater_policy']=dict(mode='downstream_free_drainage',max_generated_discharge_m3_s=1e-15,provenance='synthetic')
        with self.assertRaisesRegex(ValueError,'capacity'):forecast(c)

    def test_bypass_needs_its_own_hydraulic_input_for_correlation(self):
        c=warm_case();c['upgrades']['bypass_heat_transfer']={'method':'gnielinski_smooth'}
        with self.assertRaisesRegex(ValueError,'own explicit'):forecast(c)

    def test_local_correlation_and_declared_fixed_bypass_are_explicit(self):
        c=warm_case();c['upgrades']['toe_heat_transfer']={'method':'gnielinski_smooth'}
        r=forecast(c)
        self.assertTrue(all(x['toe_heat_transfer']['method']=='gnielinski_smooth' for x in r['series']))
        self.assertTrue(all(x['bypass_heat_transfer'] is None for x in r['series']))
        self.assertLess(abs(r['checks']['global_energy_residual_J']),1e-5)

    def test_thermal_and_geometry_time_refinement(self):
        c=warm_case();c['step_s']=60.;a=forecast(c)['totals']['face_melt_m3']
        c['step_s']=30.;b=forecast(c)['totals']['face_melt_m3']
        c['step_s']=15.;d=forecast(c)['totals']['face_melt_m3']
        self.assertLess(abs(b-d),max(abs(a-b),1e-14))

    def test_prescribed_all_cover_source_is_counted_once(self):
        c=warm_case();c['upgrades']['covered_heat_policy']='prescribed_all_cover'
        r=forecast(c)
        self.assertLess(abs(r['checks']['global_energy_residual_J']),1e-5)
        self.assertEqual(r['series'][0]['covered_heat_policy'],'prescribed_all_cover')

    def test_profile_uses_real_volume_not_uniform_routing_proxy(self):
        c=warm_case();c['physics']['face_cover_fraction']=0.
        c['periods'][0].update(duration_s=7200.,water_superheat_K=1.,face_atmospheric_net_W_m2=10.)
        c['step_s']=120.;c['upgrades'].update(face_profile_bands=32,ledge_heat_policy='adiabatic')
        r=forecast(c);last=r['series'][-1]
        self.assertTrue(last['geometry_is_routing_proxy'])
        self.assertGreater(last['undercut_diagnostic']['maximum_horizontal_undercut_m'],0.)
        self.assertNotEqual(last['profile_geometry_ledger']['physical_minus_proxy_remaining_volume_m3'],0.)
        self.assertLess(abs(r['checks']['geometry_phase_volume_residual_m3']),1e-12)
        self.assertLess(abs(r['checks']['global_energy_residual_J']),1e-5)

    def test_profile_refuses_hidden_cold_or_ledge_physics(self):
        c=warm_case();c['upgrades']['face_profile_bands']=8
        with self.assertRaisesRegex(ValueError,'ledge'):forecast(c)
        c['upgrades']['ledge_heat_policy']='adiabatic';c['periods'][0]['face_atmospheric_net_W_m2']=-1.
        with self.assertRaisesRegex(ValueError,'nonnegative'):forecast(c)

    def test_spatial_profile_and_two_column_thermal_are_coupled(self):
        c=warm_case();c['upgrades'].update(face_profile_bands=32,ledge_heat_policy='adiabatic')
        r=forecast(c)
        self.assertTrue(all(row['wet_debris_thermal']['partition'] is not None for row in r['series']))
        self.assertLess(abs(r['checks']['geometry_phase_volume_residual_m3']),1e-12)
        self.assertLess(abs(r['checks']['global_energy_residual_J']),1e-5)

    def test_initial_nonplanar_geometry_requires_an_explicit_reference_frame(self):
        c=warm_case();g=c['geometry'];mid=(g['floor_z_m']+g['crest_z_m'])/2
        c['upgrades'].update(face_profile_bands=2,ledge_heat_policy='adiabatic',initial_face_bands=[
            dict(z_low_m=g['floor_z_m'],z_high_m=mid,front_x_m=g['face_x_m']),
            dict(z_low_m=mid,z_high_m=g['crest_z_m'],front_x_m=-.1)])
        with self.assertRaisesRegex(ValueError,'provenance'):forecast(c)
        c['upgrades']['profile_provenance']='SYNTHETIC local x into ice, z up; not a field survey'
        r=forecast(c)
        self.assertGreater(r['series'][0]['undercut_diagnostic']['maximum_horizontal_undercut_m'],.09)
        self.assertLess(abs(r['checks']['geometry_phase_volume_residual_m3']),1e-12)

    def test_fixed_covered_footprint_does_not_drift_with_floor_incision(self):
        c=warm_case();g=c['geometry'];h=c['physics']['face_cover_fraction']*(g['crest_z_m']-g['floor_z_m'])
        c['upgrades'].update(fixed_cover_interval_m=[g['floor_z_m'],g['floor_z_m']+h],
            cover_interval_provenance='SYNTHETIC fixed-z patch, identical projected footprint in every comparator')
        c['periods'][0]['duration_s']=3600.;r=forecast(c)
        expected=h*g['extrusion_width_m']
        self.assertGreater(r['totals']['floor_melt_m3'],0.)
        self.assertTrue(all(abs(row['port_states']['face_debris']['area_m2']-expected)<1e-14 for row in r['series']))
        self.assertTrue(all(row['covered_z_interval_m']==tuple(c['upgrades']['fixed_cover_interval_m']) for row in r['series']))
        self.assertLess(abs(r['checks']['global_energy_residual_J']),1e-5)

    def test_package_import_isolated_from_standalone_module_cache(self):
        import importlib,importlib.util,sys
        root=Path(__file__).parent
        spec=importlib.util.spec_from_file_location('isolated_cliff_physics',root/'__init__.py',submodule_search_locations=[str(root)])
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
        library=importlib.import_module('isolated_cliff_physics.model')
        c=warm_case()
        self.assertEqual(library.forecast(c)['totals'],forecast(c)['totals'])


if __name__=='__main__':unittest.main()
