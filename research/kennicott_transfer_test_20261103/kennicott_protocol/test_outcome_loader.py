import csv,io,json,tempfile,unittest,hashlib,zlib,struct
from pathlib import Path
from unittest.mock import patch
from . import outcome_loader as m
FIELDS=['site','date_start','date_end','period_length','mass_change_kg_m2','debris_thick_mm','surface_type_start','surface_type_end','seasonal_snow_depth_start_mm','seasonal_snow_depth_end_mm','summer_snow_depth_end_mm']
def row(b,mass):
 s=lambda d:__import__('datetime').date.fromisoformat(d).strftime('%m/%d/%Y')
 return [b['site'],s(b['start']),s(b['end']),str(b['period_days']),mass,str(b['debris_mm']),'Debris','Debris','0','0','0']
def csvbytes(rows):
 f=io.StringIO();w=csv.writer(f);w.writerow(FIELDS);w.writerows(rows);return f.getvalue().encode()
def fixture(temp,role='train'):
 p=Path(temp);spec=p/'spec.json';spec.write_text(json.dumps({'native_outcome_release':m.release_plan()}));pred=p/'pred.json';pred.write_text(json.dumps({'prediction_stage':'test_fixture_no_native_values'}));receipt=p/'receipt.json'
 r={'schema':'kennicott-outcome-release-v1','access_role':role,'freeze_published':True,'predictions_published':True,'publication_reference':'synthetic-test-publication','prediction_stage':{'train':'normalized_training','evaluation':'final_evaluation'}[role],'frozen_spec_path':'spec.json','frozen_spec_sha256':m.sha(spec),'predictions_path':'pred.json','predictions_sha256':m.sha(pred),'approved_operator':m.OPERATOR};receipt.write_text(json.dumps(r));return receipt,r
class Tests(unittest.TestCase):
 def test_train_never_converts_evaluation_values(self):
  raw=csvbytes([row(m.BINDINGS['train'][0],'-100')]+[row(b,'FORBIDDEN_NONNUMERIC') for b in m.BINDINGS['evaluation']])
  x=m._select_rows(raw,m.BINDINGS['train']);self.assertEqual(x[0]['observed_net_ablation_mmwe'],100);self.assertFalse(x[0]['density_multiplier_applied'])
 def test_evaluation_never_converts_train_value(self):
  raw=csvbytes([row(m.BINDINGS['train'][0],'FORBIDDEN_NONNUMERIC')]+[row(b,str(-10-i)) for i,b in enumerate(m.BINDINGS['evaluation'])])
  x=m._select_rows(raw,m.BINDINGS['evaluation']);self.assertEqual([v['observed_net_ablation_mmwe'] for v in x],[10,11,12])
 def test_negative_target_retained(self):
  x=m._select_rows(csvbytes([row(m.BINDINGS['train'][0],'7')]),m.BINDINGS['train'])[0];self.assertEqual(x['observed_net_ablation_mmwe'],-7);self.assertTrue(x['negative_net_ablation_retained'])
 def test_duplicate_refused(self):
  r=row(m.BINDINGS['train'][0],'-100')
  with self.assertRaises(ValueError):m._select_rows(csvbytes([r,r]),m.BINDINGS['train'])
 def test_wrong_snow_surface_geometry_refused(self):
  for col,value in [('seasonal_snow_depth_end_mm','1'),('surface_type_start','Firn'),('debris_thick_mm','211'),('mass_change_kg_m2','nan')]:
   r=row(m.BINDINGS['train'][0],'-100');r[FIELDS.index(col)]=value
   with self.assertRaises(ValueError):m._select_rows(csvbytes([r]),m.BINDINGS['train'])
 def test_prepublication_gate_before_fetch(self):
  with tempfile.TemporaryDirectory() as d:
   p,r=fixture(d);r['freeze_published']=False;p.write_text(json.dumps(r))
   with patch.object(m,'_fetch_verified_native') as f:
    with self.assertRaises(ValueError):m.load_after_freeze(p,expected_role='train')
    f.assert_not_called()
 def test_role_and_stage_gate(self):
  with tempfile.TemporaryDirectory() as d:
   p,r=fixture(d,'train')
   with self.assertRaises(ValueError):m.verify_receipt(p,expected_role='evaluation')
   r['prediction_stage']='final_evaluation';p.write_text(json.dumps(r))
   with self.assertRaises(ValueError):m.verify_receipt(p,expected_role='train')
 def test_tampered_predictions_and_spec_refused(self):
  with tempfile.TemporaryDirectory() as d:
   p,r=fixture(d);(Path(d)/'pred.json').write_text('changed')
   with self.assertRaises(ValueError):m.verify_receipt(p,expected_role='train')
   p,r=fixture(d);sp=Path(d)/'spec.json';spec=json.loads(sp.read_text());spec['native_outcome_release']['unit']='m_ice';sp.write_text(json.dumps(spec));r['frozen_spec_sha256']=m.sha(sp);p.write_text(json.dumps(r))
   with self.assertRaises(ValueError):m.verify_receipt(p,expected_role='train')
 def test_mocked_role_release_no_native_access(self):
  with tempfile.TemporaryDirectory() as d:
   p,r=fixture(d);raw=csvbytes([row(m.BINDINGS['train'][0],'-123')])
   with patch.object(m,'_fetch_verified_native',return_value=raw):x=m.load_after_freeze(p,expected_role='train')
   self.assertEqual(x['rows'][0]['observed_net_ablation_mmwe'],123)
 def test_verified_member_hash_and_crc(self):
  raw=b'fixture-only';co=zlib.compressobj(wbits=-15);comp=co.compress(raw)+co.flush();name=m.MEMBER.encode();h=struct.pack('<4s5H3L2H',b'PK\x03\x04',20,0,8,0,0,zlib.crc32(raw),len(comp),len(raw),len(name),0)
  with patch.object(m,'_range',side_effect=[h,name,comp]),patch.object(m,'SOURCE_BYTES',len(raw)),patch.object(m,'SOURCE_CRC32',zlib.crc32(raw)),patch.object(m,'SOURCE_SHA256',hashlib.sha256(raw).hexdigest()):self.assertEqual(m._fetch_verified_native(),raw)
  with patch.object(m,'_range',side_effect=[h,name,comp]):
   with self.assertRaises(ValueError):m._fetch_verified_native()
if __name__=='__main__':unittest.main()
