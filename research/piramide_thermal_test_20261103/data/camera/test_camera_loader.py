import unittest,tempfile,json
from pathlib import Path
import pandas as pd
from camera_loader import score_interval_mmwe,verify_receipt,sha,OPERATOR,authorized_labels
class Checks(unittest.TestCase):
 def frame(self,a,b):return pd.DataFrame({'distance_m_ice':[a,b]},index=pd.to_datetime(['2022-12-01T12:00','2022-12-15T12:00']))
 def test_units(self):self.assertAlmostEqual(score_interval_mmwe(self.frame(1.,1.1),'2022-12-01T12:00','2022-12-15T12:00')['observed_ice_ablation_mmwe'],91.7)
 def test_negative_not_clipped(self):
  x=score_interval_mmwe(self.frame(1.1,1),'2022-12-01T12:00','2022-12-15T12:00');self.assertLess(x['observed_ice_ablation_mmwe'],0);self.assertTrue(x['negative_endpoint_difference'])
 def test_absent_endpoint(self):
  with self.assertRaises(ValueError):score_interval_mmwe(self.frame(1,2),'2022-12-01T12:30','2022-12-15T12:00')
 def test_missing_outcome(self):
  with self.assertRaises(ValueError):score_interval_mmwe(self.frame(1,float('nan')),'2022-12-01T12:00','2022-12-15T12:00')
 def test_refuses_without_freeze(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'receipt.json';p.write_text('{}')
   with self.assertRaisesRegex(ValueError,'freeze'):verify_receipt(p,d)
 def test_train_excludes_evaluation(self):
  r={'access_role':'train','intervals':[{'role':'train','tower':1,'start':'2022-12-02T03:00','end':'2022-12-16T03:00'},{'role':'evaluation','tower':1,'start':'2022-12-16T03:00','end':'2022-12-30T03:00'}]}
  labels=authorized_labels(r,1);self.assertEqual(labels,[pd.Timestamp('2022-12-02T03:00'),pd.Timestamp('2022-12-16T03:00')]);self.assertNotIn(pd.Timestamp('2022-12-30T03:00'),labels)
  with self.assertRaises(ValueError):authorized_labels(r,2)
 def test_receipt_binds_artifacts(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'spec').write_text('SYNTHETIC');(p/'predictions').write_text('SYNTHETIC');(p/'SEALED_SOURCE_MANIFEST.json').write_text('SYNTHETIC')
   r={'schema':'piramide-camera-unseal-v1','freeze_published':True,'publication_reference':'SYNTHETIC_TEST_ONLY','frozen_spec_path':str(p/'spec'),'frozen_spec_sha256':sha(p/'spec'),'predictions_path':str(p/'predictions'),'predictions_sha256':sha(p/'predictions'),'camera_manifest_sha256':sha(p/'SEALED_SOURCE_MANIFEST.json'),'approved_operator':OPERATOR,'clock_interpretation':'same_source_naive_local_no_UTC_assumption','access_role':'train','prediction_stage':'normalized_training','intervals':[{'role':'train','tower':1,'start':'2022-12-01T12:00','end':'2022-12-15T12:00'}]}
   camera_keys=['approved_operator','clock_interpretation','camera_manifest_sha256','intervals']
   (p/'spec').write_text(json.dumps({'camera_evaluation':{k:r[k] for k in camera_keys}}));r['frozen_spec_sha256']=sha(p/'spec')
   (p/'r.json').write_text(json.dumps(r));self.assertEqual(verify_receipt(p/'r.json',p)['approved_operator'],OPERATOR)
   r['access_role']='evaluation';(p/'r.json').write_text(json.dumps(r))
   with self.assertRaisesRegex(ValueError,'stage'):verify_receipt(p/'r.json',p)
   r['access_role']='train';(p/'r.json').write_text(json.dumps(r))
   (p/'predictions').write_text('TAMPERED')
   with self.assertRaisesRegex(ValueError,'hash'):verify_receipt(p/'r.json',p)
if __name__=='__main__':unittest.main()
