import csv,datetime as dt,tempfile,unittest
from pathlib import Path
from run_melt_comparison import positive_linear_integral,feature_dose

class FeatureTests(unittest.TestCase):
    def test_positive_linear_integral_crossing(self):
        self.assertEqual(positive_linear_integral(-2,2,600),300.)
        self.assertEqual(positive_linear_integral(2,-2,600),300.)
        self.assertEqual(positive_linear_integral(-2,-1,600),0.)
        self.assertEqual(positive_linear_integral(1,3,600),1200.)
    def source(self,directory):
        p=Path(directory)/'test.csv'
        p.write_text('SOURCE_LOCAL_LABEL,T\n2022-12-02T03:00:00,-2\n2022-12-02T03:10:00,2\n2022-12-02T03:20:00,4\n')
        return p
    def test_exact_clipped_support(self):
        with tempfile.TemporaryDirectory() as d:
            r=feature_dose(self.source(d),'T','2022-12-02T03:05:00','2022-12-02T03:15:00',hold='linear',max_gap_s=600)
            self.assertAlmostEqual(r['dose'],1050.)
    def test_gap_and_support_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=self.source(d)
            for start,end,gap in [('2022-12-02T02:00:00','2022-12-02T03:10:00',600),('2022-12-02T03:00:00','2022-12-02T03:20:00',300)]:
                with self.assertRaises(ValueError):feature_dose(p,'T',start,end,hold='left',max_gap_s=gap)
if __name__=='__main__':unittest.main()
