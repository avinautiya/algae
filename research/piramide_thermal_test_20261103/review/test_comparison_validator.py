import copy,itertools,math,tempfile,unittest
from pathlib import Path
from comparison_validator import validate_protocol,conservative_block_test,sign_tail,holm,write_freeze

def fixture():
 return dict(version=1,claim_scope='cross_glacier_generalization',candidate=dict(id='thermal',selected_using='development_only',all_considered_candidates=['thermal','alternative'],selection_rule='development blockMAE plus physical admissibility',code_sha256='a'*64),
  splits=dict(development_groups=['old'],evaluation_groups=['new1','new2']),targets_inspected_before_freeze=['old'],prediction_freeze_before_target_scoring=True,
  parameters=[dict(id='C',all_transitive_melt_dependencies=['old'],dependency_audit_complete=True)],datasets=[dict(version='v1',sha256='b'*64,source_url='https://example.org/source')],
  comparison=dict(target='native_ice_ablation_mmwe',same_spatial_footprint=True,same_time_support=True,same_target_quantity=True,observation_operator_verified=True,measurement_uncertainty_declared=True,all_prespecified_baselines_reported=True,baselines=['rate','PDD','EB'],primary_loss='block_mean_absolute_error',multiple_comparison_rule='holm_one_sided_sign',independence_blocks={'new1':'glacier1','new2':'glacier2'},block_definition='glacier',independence_justification_reviewed=True,minimum_practical_improvement_mmwe=10.,same_information_budget=True,baseline_sources={b:dict(source='fixed source',code_or_output_sha256='c'*64,native_ice_component_verified=True) for b in ['rate','PDD','EB']}))
class TestComparison(unittest.TestCase):
 def test_accepts_full_declaration_and_no_overwrite(self):
  p=fixture();self.assertTrue(validate_protocol(p)['ready'])
  with tempfile.TemporaryDirectory() as d:
   dest=Path(d)/'freeze.json';self.assertEqual(len(write_freeze(p,dest)),64)
   with self.assertRaises(ValueError):write_freeze(p,dest)
 def test_target_inspection_and_transitive_pooling_refused(self):
  p=fixture();p['targets_inspected_before_freeze'].append('new1');self.assertFalse(validate_protocol(p)['ready'])
  p=fixture();p['parameters'][0]['all_transitive_melt_dependencies'].append('new2');self.assertFalse(validate_protocol(p)['ready'])
 def test_different_information_or_unknown_output_component_refused(self):
  for key in ['same_information_budget','same_time_support','same_spatial_footprint']:
   p=fixture();p['comparison'][key]=False;self.assertFalse(validate_protocol(p)['ready'])
  p=fixture();p['comparison']['baseline_sources']['EB']['native_ice_component_verified']=False;self.assertFalse(validate_protocol(p)['ready'])
 def test_one_glacier_not_multiglacier(self):
  p=fixture();p['comparison']['independence_blocks']['new2']='glacier1';self.assertFalse(validate_protocol(p)['ready'])
 def test_exact_sign_tail_against_full_binary_enumeration(self):
  for n in range(1,9):
   outcomes=list(itertools.product([0,1],repeat=n))
   for k in range(n+1):self.assertEqual(sign_tail(k,n),sum(sum(o)>=k for o in outcomes)/len(outcomes))
 def test_uncertain_blocks_count_as_nonwins(self):
  rows=[dict(block=str(i),candidate_minus_baseline_loss_bounds_mmwe={'b':[-2,-1] if i<5 else [-1,1]}) for i in range(6)]
  r=conservative_block_test(rows,['b']);self.assertEqual(r['certified_block_wins']['b'],5);self.assertAlmostEqual(r['tests']['b']['raw_p'],7/64)
 def test_family_penalty_and_pseudoreplication(self):
  def rows(n):return [dict(block=str(i),candidate_minus_baseline_loss_bounds_mmwe={b:[-2,-1] for b in ['x','y','z']}) for i in range(n)]
  self.assertFalse(conservative_block_test(rows(5),['x','y','z'])['all_baselines_pass'])
  self.assertTrue(conservative_block_test(rows(6),['x','y','z'])['all_baselines_pass'])
  duplicated=rows(6);duplicated[-1]['block']='0'
  with self.assertRaises(ValueError):conservative_block_test(duplicated,['x','y','z'])
if __name__=='__main__':unittest.main()
