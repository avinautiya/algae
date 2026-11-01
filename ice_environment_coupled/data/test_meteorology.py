import tempfile,unittest
from pathlib import Path
from audit_meteorology import FIELDS,load

class MeteorologyTests(unittest.TestCase):
    def write(self,root,row):
        p=Path(root)/'test.tsv';p.write_text('\t'.join(FIELDS)+'\n'+'\t'.join(map(str,row))+'\n');return p
    def test_missing_and_negative_are_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            row=[2013,128,100,128+1/24,1,0,3,80,-2,0,270,0,'NaN',0,0]
            r=load(self.write(d,row))[0]
            self.assertTrue(r['negative_shortwave']);self.assertIsNone(r['values']['LWOUT'])
            self.assertEqual(r['timestamp_local_label'],'2013-05-08T01:00:00')
    def test_bad_timestamp_is_not_guessed(self):
        with tempfile.TemporaryDirectory() as d:
            row=[2013,128,2460,128,1,0,3,80,0,0,270,0,0,0,0]
            with self.assertRaises(ValueError):load(self.write(d,row))
if __name__=='__main__':unittest.main()
