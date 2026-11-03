import copy,unittest
from score_trace import score
SPEC=dict(pre_start_s=0.,pulse_start_s=2.,pulse_end_s=4.,post_end_s=8.,contact_threshold=.5,maximum_sampling_gap_s=2.)
def trace(vals):return [dict(time_s=float(i),toe_fraction=v,delivered_heat_W=2.*i) for i,v in enumerate(vals)]
class Tests(unittest.TestCase):
 def test_maintained_not_reconnected(self):
  r=score(trace([1.]*9),SPEC);self.assertEqual(r['classification'],'maintained_contact');self.assertFalse(r['reconnected_after_loss']);self.assertEqual(r['postpulse']['signed_delivered_heat_J'],48.)
 def test_loss_no_reconnection(self):self.assertEqual(score(trace([1.,1.,0.,0.,0.,0.,0.,0.,0.]),SPEC)['classification'],'lost_without_reconnection')
 def test_sustained_recovery(self):self.assertEqual(score(trace([1.,0.,0.,0.,1.,1.,1.,1.,1.]),SPEC)['classification'],'sustained_postpulse_recovery')
 def test_reloss(self):self.assertEqual(score(trace([1.,0.,0.,1.,1.,0.,0.,0.,0.]),SPEC)['classification'],'reconnected_then_relost')
 def test_late_recovery_not_full_window(self):self.assertEqual(score(trace([1.,0.,0.,0.,0.,0.,1.,1.,1.]),SPEC)['classification'],'temporary_or_incomplete_reconnection')
 def test_initial_dry_not_genuine_loss(self):
  r=score(trace([0.,0.,0.,0.,1.,1.,1.,1.,1.]),SPEC);self.assertFalse(r['genuine_loss']);self.assertEqual(r['classification'],'initially_disconnected_recovered')
 def test_fractional_crossing_duration(self):
  r=score(trace([1.,1.,1.,1.,0.,1.,1.,1.,1.]),SPEC);self.assertAlmostEqual(r['postpulse']['dry_seconds'],.5)
 def test_identical_history_grouped(self):
  a=trace([1.,1.,0.,0.,0.,0.,0.,0.,0.]);b=trace([1.,1.,0.,0.,1.,1.,1.,1.,1.]);self.assertEqual(score(a,SPEC)['sampled_preloss_history_sha256'],score(b,SPEC)['sampled_preloss_history_sha256'])
 def test_malformed_refused(self):
  for bad in ['gap','nan','time','fraction','bool']:
   x=trace([1.]*9)
   if bad=='gap':x.pop(4);x.pop(4)
   elif bad=='nan':x[1]['delivered_heat_W']=float('nan')
   elif bad=='time':x[1]['time_s']=0.
   elif bad=='fraction':x[1]['toe_fraction']=1.1
   else:x[1]['toe_fraction']=True
   with self.assertRaises(ValueError):score(x,SPEC)
 def test_signed_heat_not_clipped(self):
  x=trace([1.]*9)
  for r in x:r['delivered_heat_W']=-1.
  self.assertEqual(score(x,SPEC)['postpulse']['signed_delivered_heat_J'],-4.)
if __name__=='__main__':unittest.main()
