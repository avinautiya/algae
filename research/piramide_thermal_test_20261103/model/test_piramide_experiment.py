import copy,csv,datetime as dt,hashlib,json,tempfile,unittest
from pathlib import Path
from piramide_experiment import prepare,predict_frozen,check_config


ROOT=Path(__file__).resolve().parent


class FrozenColumnCampaignTests(unittest.TestCase):
    def fixture(self,directory):
        p=Path(directory);a=dt.datetime(2022,12,2,3)
        times=[a+dt.timedelta(seconds=1800*i) for i in range(9)]
        source=p/'process.csv'
        with source.open('w') as f:
            writer=csv.writer(f);writer.writerow(['SOURCE_LOCAL_LABEL','DT1','DT2','DT3','DT4','DT5'])
            for t in times:writer.writerow([t.isoformat()]+[5*(.52-z)/.5 for z in [.02,.125,.25,.375,.5]])
        h=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
        manifest=p/'process_manifest.json';manifest.write_text(json.dumps(dict(checks=[dict(family='debtemp',tower='Tower1',rows=9,sha256=h(source))])))
        clocks=p/'camera_times.json';clocks.write_text(json.dumps([t.isoformat() for t in times]))
        camera=p/'camera_manifest.json';camera.write_text('{"sealed":true}')
        c=json.loads((ROOT/'PIRAMIDE_EXPERIMENT_SPEC.json').read_text())
        c.update(process_csv=str(source),process_manifest=str(manifest),camera_times_only=str(clocks),camera_metadata=str(camera),
                 start=times[0].isoformat(),split=times[4].isoformat(),end=times[-1].isoformat(),
                 alpha_grid_m2_s=[2e-7,5e-7],cells=8,dt_seconds=300,geometries=[c['geometries'][1]])
        c['camera_evaluation']['camera_manifest_sha256']=h(camera)
        c['camera_evaluation']['intervals']=[dict(tower=1,role='train',start=c['start'],end=c['split']),dict(tower=1,role='evaluation',start=c['split'],end=c['end'])]
        config=p/'spec.json';config.write_text(json.dumps(c));return config,c

    def test_frozen_prediction_without_camera_outcomes(self):
        with tempfile.TemporaryDirectory() as td:
            config,c=self.fixture(td);out=Path(td)/'freeze';f=prepare(config,out)
            self.assertEqual(f['camera_evaluation'],c['camera_evaluation'])
            result=predict_frozen(out,'a'*40)
            self.assertFalse(result['camera_outcomes_loaded'])
            prediction=json.loads(Path(result['output']).read_text())
            candidate=[x for x in prediction['results'] if x['method']=='training_temperature_forward_fit']
            self.assertEqual(len(candidate),2)
            for row in candidate:
                self.assertEqual(len(row['normalized_heat_trace']),8)
                self.assertEqual(len(row['predicted_internal_temperature_C']),9)
                self.assertTrue(row['no_refreezing_melt_conversion_eligible'])
                self.assertEqual(row['inference']['training_rows'],4)
                self.assertLess(row['inference']['training_end_s'],row['inference']['evaluation_start_s'])
            with self.assertRaises(ValueError):predict_frozen(out,'a'*40)

    def test_source_mutation_refuses(self):
        with tempfile.TemporaryDirectory() as td:
            config,c=self.fixture(td);out=Path(td)/'freeze';prepare(config,out)
            Path(c['process_csv']).write_text('changed')
            with self.assertRaises(ValueError):predict_frozen(out,'a'*40)

    def test_geometry_phase_clock_and_camera_endpoint_gates(self):
        with tempfile.TemporaryDirectory() as td:
            config,c=self.fixture(td)
            for kind in ['geometry','clock','camera','outcomes','primary']:
                bad=copy.deepcopy(c)
                if kind=='geometry':bad['geometries'][0]['ice_depth_m']=.8
                if kind=='clock':bad['start']+='Z'
                if kind=='camera':bad['camera_evaluation']['intervals'][0]['start']='2022-01-01T00:00:00'
                if kind=='outcomes':bad['camera_outcomes_opened_before_freeze']=True
                if kind=='primary':bad['primary_geometry']='pick_best_melt'
                with self.assertRaises(ValueError):check_config(bad)


if __name__=='__main__':unittest.main()
