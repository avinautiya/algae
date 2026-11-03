"""Synthetic fixtures only; never loads native outcomes."""
import unittest,tempfile,json,types
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from camera_loader import OPERATOR,sha
from postscore_source_audit import audit
class Checks(unittest.TestCase):
 def test_only_train_endpoints_and_units(self):
  with tempfile.TemporaryDirectory() as td:
   p=Path(td);(p/'Ablation_Tower1.mat').write_bytes(b'SYNTHETIC')
   (p/'SEALED_SOURCE_MANIFEST.json').write_text(json.dumps({'members':[{'file':'Ablation_Tower1.mat','sha256':sha(p/'Ablation_Tower1.mat')}]}))
   intervals=[{'tower':1,'role':'train','start':'2022-12-02T03:00','end':'2022-12-16T03:00'},{'tower':1,'role':'evaluation','start':'2022-12-16T03:00','end':'2022-12-30T03:00'}]
   ce={'approved_operator':OPERATOR,'clock_interpretation':'same_source_naive_local_no_UTC_assumption','camera_manifest_sha256':sha(p/'SEALED_SOURCE_MANIFEST.json'),'intervals':intervals}
   (p/'spec').write_text(json.dumps({'camera_evaluation':ce}));(p/'predictions').write_text('SYNTHETIC')
   r={'schema':'piramide-camera-unseal-v1','freeze_published':True,'publication_reference':'SYNTHETIC_TEST_ONLY','frozen_spec_path':str(p/'spec'),'frozen_spec_sha256':sha(p/'spec'),'predictions_path':str(p/'predictions'),'predictions_sha256':sha(p/'predictions'),'access_role':'train','prediction_stage':'normalized_training',**ce}
   (p/'receipt').write_text(json.dumps(r))
   dates=pd.to_datetime(['2022-12-02T03:00','2022-12-16T03:00','2022-12-30T03:00']);epoch=dates.values.astype('datetime64[ms]').astype('int64')
   rowtimes=types.SimpleNamespace(properties={'data':epoch,'tz':np.array([])})
   a=np.empty((1,1),dtype=[('varNames','O'),('rowTimes','O'),('data','O')]);names=np.empty((1,2),object);names[0,0]=np.array('distance');names[0,1]=np.array('ablation')
   cells=np.empty((1,2),object);cells[0,0]=np.array([1,1.1,99]);cells[0,1]=np.array([0,.0917,-99])
   a['varNames'][0,0]=names;a['rowTimes'][0,0]=rowtimes;a['data'][0,0]=cells
   fake={'tracks':types.SimpleNamespace(properties={'any':a})}
   with patch('matio.load_from_mat',return_value=fake):report=audit(p,1,p/'receipt')
   self.assertEqual(len(report['endpoints']),2);self.assertEqual(len(report['intervals']),1)
   x=report['intervals'][0];self.assertAlmostEqual(x['independent_ice_conversion_mmwe'],91.7);self.assertAlmostEqual(x['provided_delta_over_distance_delta'],.917)
   self.assertFalse(x['operator_changed']);self.assertNotIn('2022-12-30',json.dumps(report))
if __name__=='__main__':unittest.main()
