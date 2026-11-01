import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

try:
    from .forecast import forecast,InputError,load_document
except ImportError:
    from forecast import forecast,InputError,load_document


def fixture():
    return dict(schema_version=1,evidence_label='SYNTHETIC TEST ONLY',
        parameters=dict(area_m2=10,coefficient=.7,exponent=1.5,bypass_m=0,
                        toe_m=.2,initial_stage_m=.05,peak_m3_s=.1),
        budgets=dict(available_volume_m3=3,available_time_s=30))


class Tests(unittest.TestCase):
    def test_known_scalar_prediction(self):
        r=forecast(fixture())
        self.assertEqual(r['nominal']['classification'],'feasible_static_model')
        self.assertAlmostEqual(r['nominal']['minimum_time_s'],23.632193139860757)
        self.assertAlmostEqual(r['nominal']['minimum_volume_m3'],2.3632193139860758)

    def test_peak_and_budget_failures_are_distinct(self):
        d=fixture();d['parameters']['peak_m3_s']=.01
        p=forecast(d)['nominal']
        self.assertEqual(p['failure_reasons'],['peak_at_or_below_static_activation_discharge'])
        self.assertIsNone(p['minimum_time_s'])
        d=fixture();d['budgets']=dict(available_time_s=1,available_volume_m3=1)
        self.assertEqual(set(forecast(d)['nominal']['failure_reasons']),
                         {'insufficient_deadline','insufficient_total_inflow_volume'})

    def test_initial_contact_needs_no_extra_water(self):
        d=fixture();d['parameters']['initial_stage_m']=.2
        d['budgets']=dict(available_time_s=0,available_volume_m3=0)
        p=forecast(d)['nominal']
        self.assertEqual(p['classification'],'initially_in_contact')
        self.assertEqual(p['minimum_time_s'],0)

    def test_corner_envelope_includes_each_combination(self):
        d=fixture();d['parameter_bounds']={'area_m2':[8,12],'peak_m3_s':[.05,.15]}
        s=forecast(d)['corner_sensitivity']
        self.assertEqual(s['count'],4)
        self.assertEqual({(r['parameters']['area_m2'],r['parameters']['peak_m3_s']) for r in s['samples']},
                         {(8,.05),(8,.15),(12,.05),(12,.15)})
        self.assertEqual(s['ranges']['minimum_time_s']['nonfinite_or_unreachable_count'],2)
        self.assertIn('not a confidence',s['interpretation'])

    def test_duplicate_bound_endpoint_avoids_duplicate_corner(self):
        d=fixture();d['parameter_bounds']={'area_m2':[10,10]}
        self.assertEqual(forecast(d)['corner_sensitivity']['count'],1)

    def test_schema_rejects_bad_values(self):
        mutations=[lambda d:d.update(schema_version=True),
                   lambda d:d.update(schema_version=2),
                   lambda d:d.update(unknown=1),
                   lambda d:d['parameters'].update(area_m2=True),
                   lambda d:d['parameters'].update(area_m2=float('nan')),
                   lambda d:d['parameters'].update(coefficient=-1),
                   lambda d:d['parameters'].update(peak_m3_s='0.1'),
                   lambda d:d['budgets'].update(available_time_s=-1),
                   lambda d:d.update(parameter_bounds={'area_m2':[12,8]}),
                   lambda d:d.update(parameter_bounds={'area_m2':[1,2]}),
                   lambda d:d.update(parameter_bounds={'not_parameter':[0,1]}),
                   lambda d:d.update(evidence_label='')]
        for mutate in mutations:
            d=fixture();mutate(d)
            with self.subTest(d=d),self.assertRaises(InputError):forecast(d)

    def test_provenance_stable_order_and_changes_with_data(self):
        d=fixture();a=forecast(d)['provenance']
        b=forecast(dict(reversed(list(d.items()))))['provenance']
        self.assertEqual(a,b)
        d['parameters']['area_m2']=11
        self.assertNotEqual(a['input_canonical_sha256'],forecast(d)['provenance']['input_canonical_sha256'])
        self.assertEqual(len(a['kernel_sha256']),64)

    def test_standard_json_even_unreachable(self):
        d=fixture();d['parameters']['peak_m3_s']=.01
        text=json.dumps(forecast(d),allow_nan=False)
        self.assertNotIn('Infinity',text);self.assertNotIn('NaN',text)

    def test_loader_rejects_duplicates_and_nonstandard_constants(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'bad.json'
            for content in ('{"x":1,"x":2}','{"x":NaN}','{"x":Infinity}'):
                p.write_text(content)
                with self.assertRaises(InputError):load_document(p)

    def test_cli_good_and_bad_exit_contract(self):
        script=Path(__file__).with_name('forecast.py')
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'input.json';out=Path(temp)/'forecast.json'
            p.write_text(json.dumps(fixture()))
            result=subprocess.run([sys.executable,str(script),str(p),'--output',str(out)],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(json.loads(out.read_text())['nominal']['classification'],'feasible_static_model')
            original=p.read_bytes()
            same=subprocess.run([sys.executable,str(script),str(p),'--output',str(p)],capture_output=True,text=True)
            self.assertEqual(same.returncode,2)
            self.assertEqual(p.read_bytes(),original)
            p.write_text('{"schema_version":true}')
            result=subprocess.run([sys.executable,str(script),str(p)],capture_output=True,text=True)
            self.assertEqual(result.returncode,2)
            self.assertEqual(json.loads(result.stderr)['error'],'InputError')
            self.assertEqual(result.stdout,'')


if __name__=='__main__':unittest.main()
