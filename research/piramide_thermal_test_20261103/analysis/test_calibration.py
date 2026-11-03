import unittest
from calibration import fit_training_gain,predict_from_frozen_gain

class CalibrationTests(unittest.TestCase):
    def rows(self):return [dict(start_s=0,end_s=1,observed_mmwe=-1,dose=2),dict(start_s=1,end_s=2,observed_mmwe=7,dose=4)]
    def fit(self,rows=None):return fit_training_gain(self.rows() if rows is None else rows,train_end_s=2,evaluation_start_s=3,dose_quantity='normalized conductive heat')
    def test_signed_increment_preserved(self):
        self.assertEqual(self.fit()['gain_mmwe_per_dose'],1.)
        self.assertEqual(self.fit()['training_total_mmwe'],6.)
    def test_no_evaluation_infit(self):
        r=self.rows();r[-1]['end_s']=3
        with self.assertRaises(ValueError):self.fit(r)
    def test_outcomes_refused_inprediction(self):
        with self.assertRaises(ValueError):predict_from_frozen_gain(self.fit(),[dict(start_s=3,end_s=4,dose=2,observed_mmwe=99)],evaluation_start_s=3)
    def test_forwardprediction(self):
        out=predict_from_frozen_gain(self.fit(),[dict(start_s=3,end_s=4,dose=5)],evaluation_start_s=3)
        self.assertEqual(out[0]['predicted_mmwe'],5.)
    def test_bad_dose_notclipped(self):
        r=self.rows();r[0]['dose']=-1
        with self.assertRaises(ValueError):self.fit(r)
if __name__=='__main__':unittest.main()
