import hashlib,json,tempfile,unittest
from pathlib import Path
from preflight import verify_hashes,check_scope,run,bounded
class Tests(unittest.TestCase):
 def test_scope_guard(self):
  with self.assertRaises(ValueError):check_scope('physical_bare',True)
  self.assertFalse(check_scope('conditional_quasi',True)['environmental_validation'])
 def test_hash_and_tamper(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'a').write_bytes(b'x');(p/'DELIVERY_MANIFEST.json').write_text(json.dumps({'sha256':{'a':hashlib.sha256(b'x').hexdigest()}}));verify_hashes(p);(p/'a').write_bytes(b'y')
   with self.assertRaises(ValueError):verify_hashes(p)
 def test_manifest_traversal_refused(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'DELIVERY_MANIFEST.json').write_text(json.dumps({'sha256':{'../escape':'bad'}}))
   with self.assertRaises(ValueError):verify_hashes(p)
 def test_existing_output_refused(self):
  with tempfile.TemporaryDirectory() as d:
   with self.assertRaises(ValueError):run(Path(d)/'repo',Path(d),'physical_bare',False,1,0)
 def test_timeout(self):
  import sys
  with tempfile.TemporaryDirectory() as d:
   r=bounded([sys.executable,'-c','import time;time.sleep(5)'],d,.05,Path(d)/'log');self.assertTrue(r['timed_out'])
if __name__=='__main__':unittest.main()
