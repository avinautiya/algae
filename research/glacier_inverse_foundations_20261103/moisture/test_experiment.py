import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import pandas as pd
from experiment import prepare,normalized,calibrate,sha,MODELS
class ExperimentTests(unittest.TestCase):
    def test_source_feature_to_prediction_pipeline_without_native_camera_access(self):
        with tempfile.TemporaryDirectory()as temp:
            folder=Path(temp)/'synthetic_pipeline';prepare(folder)
            normalized(folder,'a'*40)
            receipt=folder/'TRAIN_TEST_RECEIPT.json'
            receipt.write_text(json.dumps(dict(access_role='train',predictions_sha256=sha(folder/'NORMALIZED_PREDICTIONS.json'),
                frozen_spec_sha256=sha(folder/'FREEZE.json'),approved_operator='exposed_stake_distance_endpoint_difference_m_ice_times_917_to_mmwe',
                camera_manifest_sha256='synthetic-test-only')))
            synthetic=pd.DataFrame({'distance_m_ice':[.1,.2]},index=pd.to_datetime(['2022-12-02T03:00:00','2022-12-16T03:00:00']))
            with patch('experiment.load_after_freeze',return_value=synthetic)as reader:
                result=calibrate(folder,receipt)
            self.assertEqual(reader.call_count,1)
            self.assertEqual(result['models'],9)
            p=json.loads((folder/'PREDICTIONS.json').read_text())
            self.assertEqual(set(p['models']),set(MODELS))
            self.assertAlmostEqual(p['local_train_target']['observed_ice_ablation_mmwe'],91.7)
            self.assertEqual(len(p['development_parameters']),5)
            self.assertTrue((folder/'PREDICTION_SEAL.json').exists())
            self.assertTrue(p['heldout_targets_seen_before_prediction_freeze']is False)
if __name__=='__main__':unittest.main()
