import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from spatial_forecast import ASSUMPTIONS, forecast, strict_json, volume_bound


def example():
    return dict(schema_version=1,elevation_m=[[1.,1.],[1.,1.],[1.,1.]],
                dx_m=1.,dy_m=1.,minimum_cross_section_width_m=.1,
                level_bounds_m={'lower':0.,'upper':0.},
                usable_heat_bounds_J={'lower':0.,'upper':0.},
                assumptions={k:True for k in ASSUMPTIONS})


class Tests(unittest.TestCase):
    def test_narrow_trench_is_not_full_cell(self):
        r=forecast(example())
        self.assertAlmostEqual(r['optimistic_ice_volume_bound_m3']['value'],.3)
        self.assertAlmostEqual(r['diagnostic_whole_cell_path_volume_m3']['value'],3.)
        self.assertTrue(r['reconnection_excluded_under_all_declared_bounds'])

    def test_stage_uncertainty_changes_robustness(self):
        x=example();x['level_bounds_m']['upper']=1.
        r=forecast(x)
        self.assertFalse(r['reconnection_excluded_under_all_declared_bounds'])
        self.assertEqual(r['classification'],'exclusion_depends_on_declared_values')

    def test_optimistic_elevation_and_energy(self):
        x=example();x['elevation_bounds_m']={'lower':[[0.,0.]]*3,'upper':x['elevation_m']}
        self.assertFalse(forecast(x)['reconnection_excluded_under_all_declared_bounds'])

    def test_zero_width_does_not_certify_continuum(self):
        x=example();x['minimum_cross_section_width_m']=0.
        self.assertFalse(forecast(x)['reconnection_excluded_under_all_declared_bounds'])
        x['resolved_cell_route']=True
        self.assertTrue(forecast(x)['reconnection_excluded_under_all_declared_bounds'])

    def test_walls_standard_json(self):
        x=example();x['blocked']=[[False,False],[True,True],[False,False]]
        r=forecast(x)
        self.assertTrue(r['optimistic_ice_volume_bound_m3']['infinite'])
        self.assertTrue(r['reconnection_excluded_under_all_declared_bounds'])
        self.assertNotIn('Infinity',json.dumps(r,allow_nan=False))

    def test_no_ice_above_water_and_no_heat_is_inconclusive(self):
        x=example();x['level_bounds_m']={'lower':2.,'upper':2.}
        self.assertEqual(forecast(x)['classification'],'not_excluded_for_any_declared_box_value')

    def test_equality_is_not_strict_exclusion(self):
        x=example();r=forecast(x);e=r['required_latent_energy_enclosure_J']['lower']['value']
        x['usable_heat_bounds_J']={'lower':e,'upper':e}
        self.assertFalse(forecast(x)['reconnection_excluded_under_all_declared_bounds'])

    def test_input_hash_stable_key_order(self):
        x=example();y=dict(reversed(list(x.items())))
        self.assertEqual(forecast(x)['input_sha256'],forecast(y)['input_sha256'])

    def test_invalid(self):
        bad=[]
        for value in (float('nan'),float('inf'),True,'1'):
            x=example();x['dx_m']=value;bad.append(x)
        x=example();x['blocked']=[[0,0]]*3;bad.append(x)
        x=example();x['elevation_m']=[[1],[2,3]];bad.append(x)
        x=example();x['assumptions'][ASSUMPTIONS[0]]=False;bad.append(x)
        x=example();x['level_bounds_m']={'lower':1,'upper':0};bad.append(x)
        x=example();x['elevation_bounds_m']={'lower':[[2,2]]*3,'upper':[[3,3]]*3};bad.append(x)
        x=example();x['minimum_cross_section_width_m']=-1;bad.append(x)
        for x in bad:
            with self.subTest(x=x):
                with self.assertRaises(ValueError):forecast(x)

    def test_strict_json(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'input.json'
            for content in ('{"a":NaN}','{"a":1,"a":2}'):
                p.write_text(content)
                with self.assertRaises(ValueError):strict_json(p)

    def test_monotone_enclosure_random_boxes(self):
        rng=np.random.default_rng(71)
        for _ in range(15):
            z=rng.uniform(-1,1,(4,5));lo=z-.2;hi=z+.2
            low=volume_bound(lo,.7,.3,.2,.35,np.zeros(z.shape,bool),False)['volume']
            high=volume_bound(hi,.1,.3,.2,.35,np.zeros(z.shape,bool),False)['volume']
            for _ in range(10):
                zz=rng.uniform(lo,hi);h=rng.uniform(.1,.7)
                value=volume_bound(zz,h,.3,.2,.35,np.zeros(z.shape,bool),False)['volume']
                self.assertLessEqual(low,value+1e-12);self.assertLessEqual(value,high+1e-12)

    def test_numeric_overflow_is_not_physical_exclusion(self):
        x=example();x['dx_m']=1e308
        with self.assertRaises(ValueError):forecast(x)


if __name__=='__main__':unittest.main()
