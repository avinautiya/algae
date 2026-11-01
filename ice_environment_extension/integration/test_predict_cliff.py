import copy
import json
from pathlib import Path
import unittest
from predict_cliff import predict

class PredictorTests(unittest.TestCase):
    def setUp(self):self.c=json.loads((Path(__file__).parent/'synthetic_energy_packets.json').read_text())

    def test_total_volume_from_independent_energy_sum(self):
        result=predict(self.c)['cumulative']
        self.assertAlmostEqual(result['total_melt_m3'],3e6/(917*334000))
        self.assertAlmostEqual(result['polygon_volume_loss_m3'],result['total_melt_m3'],places=12)
        self.assertGreater(result['face_melt_m3'],result['floor_melt_m3'])

    def test_no_latent_heat_no_geometric_change(self):
        c=copy.deepcopy(self.c)
        for r in c['intervals']:r['face_net_latent_J']=r['floor_net_latent_J']=0.
        result=predict(c)
        self.assertEqual(result['cumulative']['total_melt_m3'],0.)
        self.assertEqual(result['series'][-1]['geometry'],c['geometry'])

    def test_provenance_and_time_required(self):
        c=copy.deepcopy(self.c);c['energy_provenance']=''
        with self.assertRaises(ValueError):predict(c)
        c=copy.deepcopy(self.c);c['intervals'][1]['end_time_s']=100.
        with self.assertRaises(ValueError):predict(c)

if __name__=='__main__':unittest.main()
