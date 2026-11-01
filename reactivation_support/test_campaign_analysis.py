import unittest
from .campaign_analysis import summarize,contact_intervals,refinement_diagnostic

def data(t,f,removal=None):return dict(days=t[-1],removal_day=removal,records=[dict(day=x,f_toe=y,Q_toe=.07*y,Q_byp=.07*(1-y)) for x,y in zip(t,f)])

class AnalysisTests(unittest.TestCase):
    def test_initially_inactive_is_not_silently_connected(self):
        r=summarize(data([0,1,2],[0,0,0]),threshold=.1,persistence_hours=12)
        self.assertFalse(r['initially_active']);self.assertEqual(r['total_observed_contact_days'],0)
    def test_crossing_time_and_observation_bracket(self):
        r=summarize(data([0,1,2],[.2,0,0]),threshold=.1,persistence_hours=12)
        self.assertEqual(r['first_loss']['day'],.5);self.assertEqual(r['first_loss']['observation_bracket_days'],[0.,1.])
    def test_flood_only_rewetting_not_ordinary_reconnection(self):
        r=summarize(data([0,1,2,3,4,5],[.2,0,.2,0,0,0]),threshold=.1,persistence_hours=12,ordinary_from_day=3.)
        self.assertIsNone(r['sustained_reconnection_start_day'])
    def test_rewet_before_return_can_qualify_if_it_persists(self):
        r=summarize(data([0,1,2,3,4,5],[.2,0,.2,.2,.2,.2]),threshold=.1,persistence_hours=24,ordinary_from_day=3.)
        self.assertEqual(r['sustained_reconnection_start_day'],3.)
    def test_never_lost_not_labeled_rescued(self):
        r=summarize(data([0,1,2],[.2,.2,.2],removal=1),threshold=.1,persistence_hours=12)
        self.assertIsNone(r['first_loss']);self.assertIsNone(r['sustained_reconnection_start_day'])
    def test_insufficient_followup_retained(self):
        r=summarize(data([0,1,2],[.2,0,.2]),threshold=.1,persistence_hours=24,ordinary_from_day=1.)
        self.assertTrue(r['insufficient_followup']);self.assertIsNone(r['sustained_reconnection_start_day'])
    def test_duplicate_times_rejected(self):
        with self.assertRaises(ValueError):contact_intervals([0,0],[0,.2],.1)
    def test_nan_rejected(self):
        with self.assertRaises(ValueError):contact_intervals([0,1],[0,float('nan')],.1)
    def test_refinement_recovers_known_order(self):
        r=refinement_diagnostic([1.04,1.01,1.0025])
        self.assertAlmostEqual(r['observed_order'],2.,places=10);self.assertAlmostEqual(r['extrapolated_value'],1.,places=10)
    def test_censoring_not_converted_to_event_at_horizon(self):
        self.assertEqual(refinement_diagnostic([3.,None,2.])['status'],'three uncensored values required')
    def test_oscillatory_refinement_does_not_get_error_bar(self):
        self.assertNotIn('fine_grid_estimated_error',refinement_diagnostic([1.,1.1,1.05]))
    def test_sampling_interval_recorded(self):
        r=summarize(data([0,.25,.5],[.2,.2,.2]),threshold=.1,persistence_hours=6)
        self.assertEqual(r['max_sampling_interval_hours'],6.)

class ActualTimeTests(unittest.TestCase):
    def test_uses_actual_removal_timestamp(self):
        run=data([0,1,2,3,4,5],[.2,0,.2,.2,.2,.2],removal=1.2)
        run['export']=[dict(t=2*86400,volume_m3=.1)]
        r=summarize(run,threshold=.1,persistence_hours=24)
        self.assertEqual(r['ordinary_from_day'],2)
        self.assertEqual(r['requested_removal_day'],1.2)
    def test_truncated_run_not_marked_complete(self):
        run=data([0,1,2],[.2,.2,.2]);run['days']=5.
        r=summarize(run,threshold=.1,persistence_hours=24)
        self.assertFalse(r['requested_horizon_completed'])

if __name__=='__main__':unittest.main()
