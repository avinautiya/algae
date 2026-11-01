import copy,json,unittest
from pathlib import Path
from model import forecast,contact_depth

class CoupledTests(unittest.TestCase):
    def setUp(self):self.c=json.loads((Path(__file__).parent/'example.json').read_text())

    def test_coupled_energy_geometry_and_routing(self):
        r=forecast(self.c)
        self.assertLess(abs(r['checks']['global_energy_residual_J']),1e-4)
        self.assertLess(abs(r['checks']['geometry_phase_volume_residual_m3']),1e-12)
        self.assertTrue(all(abs(row['routing_residual_m3_s'])<1e-12 for row in r['series']))
        self.assertGreater(r['totals']['fresh_ice_warming_J'],0.)
        self.assertGreater(r['totals']['face_melt_m3'],0.)
        self.assertTrue(any(row['time_s']==43200 for row in r['series']))
        self.assertFalse(next(row for row in r['series'] if row['time_s']==43200)['removed'])
        self.assertTrue(next(row for row in r['series'] if row['time_s']>43200)['removed'])

    def test_zero_flux_identity(self):
        for row in self.c['periods']:
            for k in row:
                if k!='duration_s':row[k]=0.
        r=forecast(self.c)
        self.assertEqual(r['series'][-1]['geometry'],self.c['geometry'])
        self.assertEqual(r['totals']['latent_J'],0.)

    def test_timestep_refinement_of_melt(self):
        a=forecast(self.c);c=copy.deepcopy(self.c);c['step_s']/=2;b=forecast(c)
        c['step_s']/=2;d=forecast(c)
        ea=abs(a['totals']['face_melt_m3']-b['totals']['face_melt_m3'])
        eb=abs(b['totals']['face_melt_m3']-d['totals']['face_melt_m3'])
        self.assertLess(eb,ea)

    def test_scope_refusals(self):
        c=copy.deepcopy(self.c);c['physics']['face_cover_fraction']=1.
        with self.assertRaises(ValueError):forecast(c)
        c=copy.deepcopy(self.c);c['physics']['max_meltwater_fraction']=0.
        with self.assertRaises(ValueError):forecast(c)
        c=copy.deepcopy(self.c);c['physics']['ice_subcooling_K']=True
        with self.assertRaises(ValueError):forecast(c)

    def test_downstream_rating_is_not_upstream_stage(self):
        r=forecast(self.c)
        row=r['series'][0]
        h=row['contact_depth_m'];w=.5;n=.02;s=.01
        independent=w*h*(w*h/(w+2*h))**(2/3)*s**.5/n
        self.assertAlmostEqual(independent,row['toe_flow_m3_s'],places=12)
        self.assertNotEqual(row['contact_stage_m'],row['upstream_stage_m'])

    def test_tiny_positive_period_is_not_silently_skipped(self):
        self.c['periods']=self.c['periods'][:1];self.c['periods'][0]['duration_s']=1e-12
        self.c['removal_time_s']=None
        self.assertEqual(len(forecast(self.c)['series']),1)

    def test_incompatible_head_is_refused(self):
        self.c['physics']['contact_friction_slope']=10.
        with self.assertRaisesRegex(ValueError,'more head'):forecast(self.c)

    def test_bulk_cold_energy_against_constant_flux_solution(self):
        self.c['physics']['face_cover_fraction']=0.
        self.c['physics']['initial_surface_cold_J_m2']=0.
        self.c['physics']['new_surface_cold_J_m2']=0.
        self.c['periods']=[dict(duration_s=3600.,inflow_m3_s=.07,water_superheat_K=0.,
                              face_atmospheric_net_W_m2=100.,floor_atmospheric_net_W_m2=0.,bypass_atmospheric_net_W_m2=0.)]
        self.c['removal_time_s']=None
        r=forecast(self.c);g=self.c['geometry']
        area=(g['crest_z_m']-g['floor_z_m'])*g['extrusion_width_m']
        energy=100.*(area-r['series'][0]['wet_bare_face_area_m2'])*3600.
        expected=energy/(917*(334000+2100*2))
        self.assertAlmostEqual(r['totals']['face_melt_m3'],expected,places=12)

if __name__=='__main__':unittest.main()
