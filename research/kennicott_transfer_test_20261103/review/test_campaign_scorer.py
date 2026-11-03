import unittest,copy
from campaign_scorer import score,canonical_sha
class Scorer(unittest.TestCase):
 def fixture(self):
  op={'kind':'campaign_total','quantity':'ice_ablation','unit':'cm_ice','source_operator_reference':'synthetic'}
  t={'target_operator':op,'rows':[{'case_id':'a','campaign_id':'c1','source_stake_id':'s','source_thermal_binding_id':'t1','period_start':'2014-01-01','period_end':'2014-01-03','observed':3.}]}
  p={'target_operator':copy.deepcopy(op),'models':{'m':[{'case_id':'a','source_thermal_binding_id':'t1','predicted':5.}]}}
  return p,t,{'access_role':'evaluation','source_binding_resolved':True,'frozen_predictions_sha256':canonical_sha(p),'released_typed_targets_sha256':canonical_sha(t)}
 def test_same_operator(self):
  p,t,r=self.fixture();self.assertEqual(score(p,t,r)['models']['m']['MAE_native_unit'],2)
 def test_refuse_hash(self):
  p,t,r=self.fixture();p['models']['m'][0]['predicted']=3
  with self.assertRaises(ValueError):score(p,t,r)
 def test_refuse_units(self):
  p,t,r=self.fixture();p['target_operator']['unit']='mmwe';r['frozen_predictions_sha256']=canonical_sha(p)
  with self.assertRaises(ValueError):score(p,t,r)
 def test_refuse_footprint(self):
  p,t,r=self.fixture();p['models']['m'][0]['source_thermal_binding_id']='other';r['frozen_predictions_sha256']=canonical_sha(p)
  with self.assertRaises(ValueError):score(p,t,r)
 def test_refuse_rate_inference(self):
  p,t,r=self.fixture();p['target_operator']['kind']=t['target_operator']['kind']='campaign_mean_daily_rate';r.update(frozen_predictions_sha256=canonical_sha(p),released_typed_targets_sha256=canonical_sha(t))
  with self.assertRaises(ValueError):score(p,t,r)
 def test_refuse_unbound(self):
  p,t,r=self.fixture();r['source_binding_resolved']=False
  with self.assertRaises(ValueError):score(p,t,r)
if __name__=='__main__':unittest.main()

class NetLossScorer(unittest.TestCase):
 def test_signed_net_loss_not_clipped(self):
  p,t,r=Scorer().fixture();p['target_operator']['quantity']=t['target_operator']['quantity']='net_ice_ablation_loss_conditional';t['rows'][0]['observed']=-1.;r.update(frozen_predictions_sha256=canonical_sha(p),released_typed_targets_sha256=canonical_sha(t));self.assertEqual(score(p,t,r)['models']['m']['all_rows'][0]['observed'],-1.)
