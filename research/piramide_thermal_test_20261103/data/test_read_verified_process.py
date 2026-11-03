import unittest,json,tempfile,shutil
from pathlib import Path
from read_verified_process import read_process
ROOT=Path(__file__).parent/'decoded'
class Checks(unittest.TestCase):
    def test_source_thermal_shapes(self):
        for t,n in [('Tower1',4760),('Tower2',4719),('Tower3',4907)]:
            f=read_process(ROOT,tower=t);self.assertEqual(f.shape,(n,5));self.assertEqual(list(f.columns),['DT1','DT2','DT3','DT4','DT5'])
    def test_original_labels(self):
        f=read_process(ROOT);self.assertEqual(str(f.index[0]),'2022-11-24 08:00:00');self.assertIsNone(f.index.tz)
    def test_fail_modified_csv(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for n in ['DECODE_VERIFICATION.json','debtemp_Tower1.csv']:shutil.copy2(ROOT/n,p/n)
            with (p/'debtemp_Tower1.csv').open('a') as out:out.write('\n')
            with self.assertRaisesRegex(ValueError,'hash'):read_process(p)
    def test_not_exact_temperate_basal(self):
        f=read_process(ROOT);self.assertGreater(f.DT5.max(),4);self.assertGreater((f.DT5>.4).sum(),2000)
if __name__=='__main__':unittest.main()
