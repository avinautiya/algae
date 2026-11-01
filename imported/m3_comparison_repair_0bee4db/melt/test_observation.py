import copy
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from .observation import rectangle_overlap,geometry_melt,latent_energy,volume_from_latent,profile_melt,reduced_axial_roi,export_sim_interval
from .repair_eval import repair,validate_output_path


class Tests(unittest.TestCase):
    def test_fractional_overlap_and_nonuniform_grid(self):
        weights=rectangle_overlap([0,1,3],[0,2,3],[.5,2.,1.,2.5])
        np.testing.assert_allclose(weights,[[.5,.25],[1.,.5]])
        self.assertAlmostEqual(weights.sum(),2.25)

    def test_geometry_volume_and_mask(self):
        a=np.zeros((2,2));b=a-.1
        r=geometry_melt(a,b,[0,1,2],[0,1,2],[.5,1.5,0,2],np.array([[True,False],[True,True]]))
        self.assertAlmostEqual(r['volume_m3'],.15)
        self.assertAlmostEqual(r['observed_ice_area_m2'],1.5)
        with self.assertRaises(ValueError):geometry_melt(a,b+.2,[0,1,2],[0,1,2],[0,2,0,2])

    def test_energy_volume_roundtrip(self):
        self.assertAlmostEqual(volume_from_latent(latent_energy(.123)),.123)
        self.assertEqual(latent_energy(1.),917*334000)

    def test_existing_profile_control_open_decomposition(self):
        first=[-.08]*80;last=(np.array(first)-np.linspace(.01,.04,80)).tolist()
        rows=[dict(day=0,toe_floor_x=first),dict(day=21,toe_floor_x=last)]
        common=dict(records=rows,dx_m=.025,toe_width_m=.5,reach_start_m=3.)
        all_=profile_melt(**common,roi_start_m=3.,roi_end_m=5.)
        ctrl=profile_melt(**common,roi_start_m=3.5,roi_end_m=4.5)
        left=profile_melt(**common,roi_start_m=3.,roi_end_m=3.5)
        right=profile_melt(**common,roi_start_m=4.5,roi_end_m=5.)
        self.assertAlmostEqual(all_['volume_m3'],ctrl['volume_m3']+left['volume_m3']+right['volume_m3'])
        self.assertAlmostEqual(ctrl['area_m2'],.5)
        self.assertAlmostEqual(ctrl['rounding_volume_bound_m3'],5e-6)

    def test_profile_fractional_axial_cells(self):
        rows=[dict(day=0,toe_floor_x=[0.,0.]),dict(day=1,toe_floor_x=[-.1,-.2])]
        r=profile_melt(rows,dx_m=1.,toe_width_m=.5,reach_start_m=0.,roi_start_m=.5,roi_end_m=1.25)
        self.assertAlmostEqual(r['volume_m3'],.05)

    def test_reduced_common_roi(self):
        self.assertAlmostEqual(reduced_axial_roi(.2,method='connectivity_reduced',layout='strip')['volume_m3'],.2)
        self.assertAlmostEqual(reduced_axial_roi(.2,method='connectivity_reduced',layout='corridor')['volume_m3'],.1)
        for method in ('fractional_cover','fixed_geometry','bare_reduced(evolving)','bare_fixed_geometry'):
            self.assertAlmostEqual(reduced_axial_roi(.2,method=method,layout='strip')['volume_m3'],.1)
        with self.assertRaises(ValueError):reduced_axial_roi(.2,method='connectivity_reduced',layout='strip',roi_start_from_reach_m=0.,roi_end_from_reach_m=2.)
        with self.assertRaises(ValueError):reduced_axial_roi(.2,method='connectivity_reduced',layout='porous_strip')

    def test_export_reports_floor_clipping_disagreement(self):
        before=np.zeros((4,4));after=before.copy();after[1:-1,1:-1]=-.01
        sim=SimpleNamespace(nx=2,ny=2,dx=1.,zb=after,ice=np.ones((4,4),bool),wall=np.zeros((4,4),bool),last_qacc=np.full((4,4),917*334000*.02))
        r=export_sim_interval(sim,before,roi_rectangle=[0,2,0,2],dt_s=1.)
        self.assertAlmostEqual(r['volume_m3'],.04)
        self.assertAlmostEqual(r['heat_implied_volume_m3'],.08)
        self.assertAlmostEqual(r['geometric_minus_heat_volume_m3'],-.04)

    def test_missing_raw_data_never_scores(self):
        spec={'evaluation_runs':[dict(name='E04_strip_none',layout='strip',dx=.025,removal_day=None,days=21.,dtm_h=6.)],'reduced_predictions':{}}
        with tempfile.TemporaryDirectory() as d:
            r=repair(d,spec)
        self.assertEqual(r['results'][0]['status'],'unavailable')
        self.assertEqual(r['results'][0]['comparisons'],[])

    def test_synthetic_raw_match_all_model_scopes(self):
        spec={'evaluation_runs':[dict(name='E04_strip_none',layout='strip',dx=.025,removal_day=None,days=21.,dtm_h=6.)],'reduced_predictions':{'E04_strip_none':[
            dict(method='connectivity_reduced',toe_melt_m3=.05),dict(method='fractional_cover',toe_melt_m3=.1)]}}
        raw=dict(status='complete',dx=.025,layout='strip',removal_day=None,days=21.,dtm=21600.,records=[dict(day=0,toe_floor_x=[-.08]*80),dict(day=21,toe_floor_x=[-.18]*80)],cum={'toe_reach':.1})
        with tempfile.TemporaryDirectory() as d:
            Path(d,'E04_strip_none_medium.json').write_text(json.dumps(raw));r=repair(d,spec)
        self.assertEqual(r['results'][0]['status'],'reconstructed')
        for c in r['results'][0]['comparisons']:self.assertAlmostEqual(c['error_m3'],0.)

    def test_wrong_removal_days_and_dtm_rejected(self):
        spec={'evaluation_runs':[dict(name='E06_strip_r8',layout='strip',dx=.025,removal_day=8.,days=21.,dtm_h=6.)]}
        raw=dict(status='complete',dx=.025,layout='strip',removal_day=8.,days=21.,dtm=21600.,records=[])
        with tempfile.TemporaryDirectory() as d:
            for field,bad in [('removal_day',12.),('days',20.),('dtm',10800.)]:
                wrong=copy.deepcopy(raw);wrong[field]=bad
                Path(d,'E06_strip_r8_medium.json').write_text(json.dumps(wrong))
                with self.assertRaisesRegex(ValueError,field):repair(d,spec)

    def test_truncated_or_missing_initial_profile_not_scored(self):
        spec={'evaluation_runs':[dict(name='E04_strip_none',layout='strip',dx=.025,removal_day=None,days=21.,dtm_h=6.)]}
        raw=dict(status='complete',dx=.025,layout='strip',removal_day=None,days=21.,dtm=21600.,
                 records=[dict(day=0,toe_floor_x=[-.08]*80),dict(day=21,toe_floor_x=[-.18]*80)])
        with tempfile.TemporaryDirectory() as d:
            for endpoint,time in [(0,1.),(1,20.)]:
                wrong=copy.deepcopy(raw);wrong['records'][endpoint]['day']=time
                Path(d,'E04_strip_none_medium.json').write_text(json.dumps(wrong));r=repair(d,spec)
                self.assertEqual(r['results'][0]['status'],'unavailable')
                self.assertEqual(r['results'][0]['comparisons'],[])

    def test_wrong_initial_floor_not_scored(self):
        spec={'evaluation_runs':[dict(name='E04_strip_none',layout='strip',dx=.025,removal_day=None,days=21.,dtm_h=6.)]}
        raw=dict(status='complete',dx=.025,layout='strip',removal_day=None,days=21.,dtm=21600.,
                 records=[dict(day=0,toe_floor_x=[-.10]*80),dict(day=21,toe_floor_x=[-.18]*80)])
        with tempfile.TemporaryDirectory() as d:
            Path(d,'E04_strip_none_medium.json').write_text(json.dumps(raw));r=repair(d,spec)
        self.assertEqual(r['results'][0]['status'],'unavailable')
        self.assertIn('initial toe profile',r['results'][0]['reason'])

    def test_output_cannot_replace_any_raw_input(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ev=root/'eval';ev.mkdir();spec=root/'spec.json'
            with self.assertRaises(ValueError):validate_output_path(ev,spec,ev/'E04_strip_none_medium.json')
            with self.assertRaises(ValueError):validate_output_path(ev,spec,ev/'new_output.json')
            with self.assertRaises(ValueError):validate_output_path(ev,spec,spec)
            validate_output_path(ev,spec,root/'repaired.json')

    def test_output_hardlink_to_raw_is_forbidden(self):
        import os
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);ev=root/'eval';ev.mkdir();raw=ev/'case.json';raw.write_text('{}')
            out=root/'output.json';os.link(raw,out)
            with self.assertRaisesRegex(ValueError,'aliases'):validate_output_path(ev,root/'spec.json',out)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):rectangle_overlap([0,0],[0,1],[0,1,0,1])
        with self.assertRaises(ValueError):latent_energy(-1)
        with self.assertRaises(ValueError):profile_melt([dict(day=0,toe_floor_x=[0])],dx_m=.1,toe_width_m=.5,reach_start_m=3,roi_start_m=3,roi_end_m=3.1)


if __name__=='__main__':unittest.main()
