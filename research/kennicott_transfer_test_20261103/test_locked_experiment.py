import importlib.util,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('locked',R/'run_locked_experiment.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
class LockedPipelineTests(unittest.TestCase):
 def test_end_to_end_with_mocked_source_only(self):
  with tempfile.TemporaryDirectory()as tmp:
   root=Path(tmp);d=root/'frozen';d.mkdir();models=root/'kennicott_model';models.mkdir()
   (d/'FROZEN_SPEC.json').write_text(json.dumps({'source_hashes':{}}))
   (models/'FEATURE_ONLY_PLAN_FINAL.json').write_text('{}')
   basis={'sites':{s:{'homogeneous_empirical_reference_dose':2.,'index_doses':{m:2. for m in mod.INDICES}}for s in ['P_AWS_1','P11','P12','P13']}}
   (models/'NOMINAL_REFERENCE_EXPOSURES_FINAL.json').write_text(json.dumps(basis))
   calls=[]
   def fake_source(receipt,expected_role):
    calls.append(expected_role)
    return {'rows':[{'site':'P_AWS_1','observed_net_ablation_mmwe':100.}]}if expected_role=='train'else {'rows':[{'site':s,'observed_net_ablation_mmwe':100.}for s in ['P11','P12','P13']]}
   with patch.object(mod,'ROOT',root),patch.object(mod,'load_after_freeze',side_effect=fake_source),patch.object(mod,'calibrate_capacity',return_value={'capacity_J_m3K':1.5e6}),patch.object(mod,'forward',return_value={'physical_melt_mmwe':100.,'relative_energy_residual':0.}):
    pred=mod.calibrate(d,'mock_train');self.assertEqual(calls,['train']);self.assertEqual(len(pred['models']),10)
    receipt=root/'mock_eval';receipt.write_text('mock')
    # Review code is a pure scoring helper, not an outcome source.
    import sys;sys.path.insert(0,str(R/'review'))
    result=mod.score(d,receipt);self.assertEqual(calls,['train','evaluation']);self.assertFalse(result['decision']['flag'])
 def test_tampered_source_rejected_before_release(self):
  with tempfile.TemporaryDirectory()as t:
   d=Path(t);(d/'FROZEN_SPEC.json').write_text(json.dumps({'source_hashes':{'run_locked_experiment.py':'wrong'}}))
   with patch.object(mod,'load_after_freeze')as source:
    with self.assertRaises(ValueError):mod.calibrate(d,'unused')
    source.assert_not_called()
if __name__=='__main__':unittest.main()
