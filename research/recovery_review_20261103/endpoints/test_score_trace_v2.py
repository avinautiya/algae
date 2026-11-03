import unittest
from score_trace_v2 import score,compare_frozen_prediction,sampled_history_digest
S=dict(pre_start_s=0.,pulse_start_s=10.,pulse_end_s=20.,post_end_s=30.,contact_threshold=.5,maximum_sampling_gap_s=10.)
def trace(vals):return [dict(time_s=t,toe_fraction=f,delivered_heat_W=2*f) for t,f in vals]
def run(vals):return score(trace(vals),S,trace_interpretation='piecewise_linear')
class Tests(unittest.TestCase):
 def test_initialdry_pulsewet_thenlost(self):
  r=run([(0,0),(10,0),(15,1),(20,0),(30,0)]);self.assertEqual(r['classification'],'initially_disconnected_temporarily_rewetted_then_dry');self.assertEqual(r['historical_v1_score']['classification'],'lost_without_reconnection');self.assertFalse(r['sustained_entire_postpulse_window']);self.assertEqual(r['first_threshold_wetting_s'],12.5)
 def test_neverwet(self):self.assertEqual(run([(0,0),(10,0),(20,0),(30,0)])['classification'],'initially_disconnected_never_wet')
 def test_initialdry_sustained_access(self):self.assertEqual(run([(0,0),(10,0),(20,1),(30,1)])['classification'],'initially_disconnected_sustained_postpulse_access')
 def test_initialdry_partial(self):self.assertEqual(run([(0,0),(10,0),(20,0),(30,1)])['classification'],'initially_disconnected_partial_postpulse_access')
 def test_maintained(self):self.assertEqual(run([(0,1),(10,1),(20,1),(30,1)])['classification'],'maintained_contact')
 def test_loss_never_reconnect(self):self.assertEqual(run([(0,1),(10,0),(20,0),(30,0)])['classification'],'lost_without_reconnection')
 def test_loss_then_sustained(self):self.assertEqual(run([(0,1),(10,0),(20,1),(30,1)])['classification'],'sustained_postpulse_recovery')
 def test_reloss(self):self.assertEqual(run([(0,1),(10,0),(20,1),(30,0)])['classification'],'reconnected_then_relost')
 def test_signed_heat_unchanged(self):
  r=run([(0,1),(10,1),(20,1),(30,1)]);self.assertEqual(r['whole_window']['signed_delivered_heat_J'],60.)
 def test_frozen_mismatch_stays_failed(self):
  r=compare_frozen_prediction(predicted_classes=['temporary_or_incomplete_reconnection','initially_disconnected'],historical_class='lost_without_reconnection');self.assertFalse(r['historical_label_match']);self.assertFalse(r['retrospective_relabeling_permitted'])
 def test_gap_and_nonfinite_refusal(self):
  for vals in [[(0,0),(30,1)],[(0,0),(10,float('nan')),(20,1),(30,1)]]:
   with self.assertRaises(ValueError):run(vals)
 def test_contract_required(self):
  with self.assertRaises(ValueError):score(trace([(0,0),(10,0),(20,0),(30,0)]),S,trace_interpretation='jump_exact')
 def test_history_digest_not_silent_null(self):
  r=run([(0,1),(10,0),(20,0),(30,0)]);self.assertEqual(len(sampled_history_digest(r['historical_v1_score'])),64)
  for v in [{},{'prefix_sha256':None},{'sampled_preloss_history_sha256':None}]:
   with self.assertRaises(ValueError):sampled_history_digest(v)
 def test_threshold_touch_is_not_positive_duration(self):
  r=run([(0,0),(10,0),(15,.5),(20,0),(30,0)]);self.assertFalse(r['positive_duration_contact_observed']);self.assertEqual(r['first_threshold_wetting_s'],15.)
if __name__=='__main__':unittest.main()
