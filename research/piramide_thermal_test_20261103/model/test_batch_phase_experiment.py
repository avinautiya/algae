import datetime as dt,unittest
from batch_phase_experiment import sum_interval,interval_indices


class FutureBatchIntervalTests(unittest.TestCase):
    def test_exact_mass_latent_and_normalization(self):
        start=dt.datetime(2022,12,30)
        labels=[start+dt.timedelta(seconds=1800*i) for i in range(5)]
        rows=[dict(melt_mmwe=v) for v in [1.,2.,3.,4.]]
        r=sum_interval(rows,labels,labels[1].isoformat(),labels[3].isoformat(),2e6,334000.)
        self.assertEqual(r['raw_melt_mmwe'],5.)
        self.assertEqual(r['source_native_intervals'],2)
        self.assertEqual(r['raw_latent_heat_J_m2'],5*334000.)
        self.assertEqual(r['normalized_basal_heat_K_m'],5*334000/2e6)

    def test_consecutive_blocks_no_duplicate_middle_interval(self):
        start=dt.datetime(2022,12,30)
        labels=[start+dt.timedelta(seconds=1800*i) for i in range(5)]
        rows=[dict(melt_mmwe=v) for v in [1.,2.,3.,4.]]
        first=sum_interval(rows,labels,labels[0].isoformat(),labels[2].isoformat(),1.,1.)
        second=sum_interval(rows,labels,labels[2].isoformat(),labels[4].isoformat(),1.,1.)
        self.assertEqual(first['raw_melt_mmwe']+second['raw_melt_mmwe'],10.)

    def test_missing_or_reversed_endpoint_refused(self):
        labels=[dt.datetime(2022,12,30),dt.datetime(2022,12,31)]
        with self.assertRaises(ValueError):interval_indices(labels,'2022-12-30T00:01:00','2022-12-31T00:00:00')
        with self.assertRaises(ValueError):interval_indices(labels,'2022-12-31T00:00:00','2022-12-30T00:00:00')


if __name__=='__main__':unittest.main()
