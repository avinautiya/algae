import copy
import json
import math
import unittest
import numpy as np
from .model import forecast,G


def fixture():
    route=dict(ice_z_m=0.,width_m=.5,ice_area_m2=1.,heat_transfer_W_m2_K=200.,
               debris_thickness_m=0.,debris_conductivity_W_m_K=1.)
    return dict(pond_area_m2=10.,pond_floor_m=-1.,initial_stage_m=.1,
                initial_temperature_K=.5,horizon_s=3600.,output_step_s=600.,
                routes=[copy.deepcopy(route),copy.deepcopy(route)],
                forcing=[dict(time_s=0.,inflow_m3_s=.07,temperature_K=.5)])


class Tests(unittest.TestCase):
    def test_independent_isothermal_drainage_solution(self):
        c=fixture();c['initial_temperature_K']=0.;c['forcing'][0].update(inflow_m3_s=0.,temperature_K=0.)
        r=forecast(c)
        coefficient=math.sqrt(G)*(2/3)**1.5
        expected=(.1**(-.5)+coefficient*c['horizon_s']/(2*c['pond_area_m2']))**(-2)
        self.assertAlmostEqual(r['series'][-1]['stage_m'],expected,places=8)
        self.assertTrue(all(r['checks'].values()))

    def test_symmetry_and_coupled_ledgers(self):
        r=forecast(fixture())
        for row in r['series']:
            self.assertAlmostEqual(row['route_flow_m3_s'][0],row['route_flow_m3_s'][1])
            self.assertAlmostEqual(row['ice_z_m'][0],row['ice_z_m'][1])
        self.assertGreater(sum(r['ledger']['ice_latent_energy_J']),0.)
        self.assertTrue(all(r['checks'].values()))
        self.assertGreater(sum(r['ledger']['downstream_meltwater_m3']),0.)

    def test_zero_transfer_has_no_melt(self):
        c=fixture()
        for route in c['routes']:route['heat_transfer_W_m2_K']=0.
        r=forecast(c)
        np.testing.assert_allclose(r['series'][-1]['ice_z_m'],[0.,0.],atol=1e-18,rtol=0.)
        np.testing.assert_allclose(r['ledger']['ice_latent_energy_J'],[0.,0.],atol=1e-12,rtol=0.)

    def test_independent_constant_volume_mixing_solution(self):
        c=fixture()
        for route in c['routes']:route['heat_transfer_W_m2_K']=0.
        q=math.sqrt(G)*(2/3)**1.5*.1**1.5
        c['forcing'][0].update(inflow_m3_s=q,temperature_K=1.)
        r=forecast(c)
        volume=11.
        expected=1.-.5*math.exp(-q*c['horizon_s']/volume)
        self.assertAlmostEqual(r['series'][-1]['reservoir_temperature_K'],expected,places=8)
        self.assertAlmostEqual(r['series'][-1]['reservoir_volume_m3'],volume,places=8)

    def test_exact_events_and_final_time(self):
        c=fixture();c['routes'][0]['debris_thickness_m']=.1
        c['removals']=[dict(time_s=777.,route=0)]
        c['forcing'].append(dict(time_s=1234.,inflow_m3_s=.02,temperature_K=.1))
        c['horizon_s']=3701.
        r=forecast(c);times=[row['time_s'] for row in r['series']]
        self.assertIn(777.,times);self.assertIn(1234.,times);self.assertEqual(times[-1],3701.)
        self.assertEqual(len(times),len(set(times)))

    def test_tolerance_refinement(self):
        c=fixture();r=forecast(c);c.update(rtol=1e-10,atol=1e-12);s=forecast(c)
        np.testing.assert_allclose(r['series'][-1]['ice_z_m'],s['series'][-1]['ice_z_m'],atol=1e-9,rtol=1e-6)

    def test_removal_flow_jump_and_continuous_storage(self):
        c=fixture();c['routes'][0]['debris_thickness_m']=.1;c['output_step_s']=777.
        c['forcing'].append(dict(time_s=777.,inflow_m3_s=.07,temperature_K=.5))
        baseline=forecast(c);c['removals']=[dict(time_s=777.,route=0)];removed=forecast(c)
        a=next(x for x in baseline['series'] if x['time_s']==777.)
        b=next(x for x in removed['series'] if x['time_s']==777.)
        self.assertAlmostEqual(a['stage_m'],b['stage_m'],places=10)
        self.assertAlmostEqual(a['reservoir_sensible_energy_J'],b['reservoir_sensible_energy_J'],places=5)
        np.testing.assert_allclose(a['ice_z_m'],b['ice_z_m'],atol=1e-12)
        self.assertGreater(b['route_flow_m3_s'][0],a['route_flow_m3_s'][0])

    def test_downstream_dilution_and_nonnegative_states(self):
        r=forecast(fixture())
        for row in r['series']:
            self.assertGreaterEqual(row['reservoir_sensible_energy_J'],-1e-6)
            for q,T in zip(row['downstream_route_flow_m3_s'],row['downstream_route_temperature_K']):
                self.assertGreaterEqual(q,0.)
                self.assertGreaterEqual(T,-1e-12)
                self.assertLessEqual(T,row['reservoir_temperature_K']+1e-12)

    def test_invalid_inputs(self):
        for key,value in [('pond_area_m2',True),('initial_temperature_K',float('nan')),('output_step_s',0.)]:
            c=fixture();c[key]=value
            with self.assertRaises(ValueError):forecast(c)
        c=fixture();c['unexpected']=1
        with self.assertRaises(ValueError):forecast(c)

    def test_hash_changes_on_removal(self):
        c=fixture();r=forecast(c);c['removals']=[dict(time_s=100,route=0)]
        self.assertNotEqual(r['input_sha256'],forecast(c)['input_sha256'])
        json.dumps(r,allow_nan=False)

    def test_control_floor_scope_exhaustion_refuses(self):
        c=fixture();c['pond_floor_m']=-.0001
        with self.assertRaisesRegex(RuntimeError,'fixed pond floor'):forecast(c)


if __name__=='__main__':unittest.main()
