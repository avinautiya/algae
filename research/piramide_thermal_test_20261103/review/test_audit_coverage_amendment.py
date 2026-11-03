import copy,itertools,json,tempfile,unittest
from pathlib import Path
from audit_coverage_amendment import nearest_finite,audit,sha
class TestCoverageAudit(unittest.TestCase):
 def test_nearest_finite_tie_earlier(self):
  rows=[dict(source_local_label=t,distance_isfinite=f) for t,f in [('2022-12-30T03:00:00',False),('2022-12-30T00:00:00',True),('2022-12-30T06:00:00',True),('2022-12-29T14:00:00',True)]]
  self.assertEqual(nearest_finite(rows,'2022-12-30T03:00:00'),'2022-12-30T00:00:00')
 def test_no_eligible_finite_refuses(self):
  with self.assertRaises(ValueError):nearest_finite([dict(source_local_label='2022-12-29T14:00:00',distance_isfinite=True)],'2022-12-30T03:00:00')
 def test_no_artifact_read_without_release(self):
  with self.assertRaisesRegex(ValueError,'announce evaluation release'):audit(*['/not/a/file']*8)
 def fixture(self,root):
  def write(name,obj):p=root/name;p.write_text(json.dumps(obj));return p
  oldend='2022-12-30T03:00:00';end='2022-12-30T00:00:00';split='2022-12-16T03:00:00'
  camera=dict(intervals=[dict(role='train',start='2022-12-02T03:00:00',end=split),dict(role='evaluation',start=split,end=oldend)],operator='synthetic')
  f0=dict(config=dict(end=oldend,unchanged_physics='fixed',camera_evaluation=camera),source_hashes={'solver':'x'},input_hashes={},analysis_source_hashes={},analysis_input_hashes={})
  oldfreeze=write('oldfreeze.json',f0)
  missing=dict(original_end=oldend,selected_end=end,missing_value_native_nan_confirmed=True,new_target_values_exposed=False,prediction_error_computed=False,nearby_timestamps_finiteness_only=[dict(source_local_label=end,distance_isfinite=True),dict(source_local_label=oldend,distance_isfinite=False)])
  missingpath=write('missing.json',missing);f1=copy.deepcopy(f0);f1['config']['end']=end;f1['config']['camera_evaluation']['intervals'][1]['end']=end;f1['coverage_amendment']=dict(original_freeze_sha256=sha(oldfreeze),missingness_audit_sha256=sha(missingpath),original_end=oldend,selected_end=end)
  newfreeze=write('newfreeze.json',f1);rows=[];names=[]
  for method,g,h,c in itertools.product(['forward','OLS'],['small','middle','large'],['left','right'],['dry','proxy','saturated']):
   name='__'.join([method,g,h,c]);names.append(name)
   row=dict(method=method,geometry={'id':g},hold=h,capacity_id=c,physicalcapacity_J_m3K=1e6,capacity_J_m3K=1e6,alpha_m2_s=1e-6,debris_conductivity_W_mK=1,alpha_refitted=False,training_normalized_basal_heat_K_m=1,raw_training_melt_mmwe=3,normalized_heat_trace=[dict(end=split),dict(end=end),dict(end=oldend)])
   rows.append(row)
  n0={'results':rows};n1=copy.deepcopy(n0)
  for row in n1['results']:row['normalized_heat_trace'].pop()
  names+=['constant','PDD1','PDD2'];p0=dict(primary_candidate=names[0],baseline_ids=names[-3:],secondary_model_ids=names[1:-3],models={n:[] for n in names},calibration={'known_train':10},calibrations={n:{'gain':2} for n in names},target_operator={'native':True},feature_doses={n:{'train':{'dose':3}} for n in ['positive_buried_T_DT1','positive_air_TA3_linear_gapfill18h']});p1=copy.deepcopy(p0)
  failure=dict(status='unscorable',evaluation_error_computed=False,score_written=False)
  prior=root/'prior';prior.mkdir();oldpred=prior/'PREDICTIONS.json';oldpred.write_text(json.dumps(p0))
  return [newfreeze,oldfreeze,write('newnorm.json',n1),write('oldnorm.json',n0),write('newpred.json',p1),oldpred,write('failure.json',failure),missingpath]
 def test_full_synthetic_prefix_gain_audit(self):
  with tempfile.TemporaryDirectory() as d:
   args=self.fixture(Path(d));r=audit(*args,evaluation_released=True);self.assertTrue(r['all_checks_pass'],[c for c in r['checks'] if not c['pass_']]);self.assertEqual(r['endpoint_shift_s'],-10800)
 def test_train_gain_change_detected(self):
  with tempfile.TemporaryDirectory() as d:
   args=self.fixture(Path(d));p=args[4];x=json.loads(p.read_text());next(iter(x['calibrations'].values()))['gain']=3;p.write_text(json.dumps(x));self.assertFalse(audit(*args,evaluation_released=True)['all_checks_pass'])
if __name__=='__main__':unittest.main()
