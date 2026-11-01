import unittest
from recovery import classify,history_key,classify_run


class Tests(unittest.TestCase):
    def test_maintained_is_not_reconnection(self):
        r=classify([0,1,2,3],[.1]*4,actual_removal_day=1)
        self.assertEqual(r['classification'],'maintained_contact')
        self.assertFalse(r['loss_before_actual_removal'])
    def test_true_reconnection(self):
        r=classify([0,1,2,3,4],[.1,0,.1,.1,.1],actual_removal_day=1.5)
        self.assertEqual(r['classification'],'genuine_reconnection')
    def test_terminal_wet_is_not_automatically_sustained(self):
        r=classify([0,1,2,2.5],[.1,0,0,.1],actual_removal_day=1.5)
        self.assertEqual(r['classification'],'unresolved_terminal_wet_episode')
    def test_persistence_is_counted_after_actual_removal(self):
        r=classify([0,1,2,3],[.1,0,.1,.1],actual_removal_day=2.5)
        self.assertNotEqual(r['classification'],'genuine_reconnection')
    def test_recovery_known_before_removal_is_not_attributed_to_it(self):
        r=classify([0,1,2,3,4],[.1,0,.1,.1,.1],actual_removal_day=3)
        self.assertEqual(r['classification'],'reconnected_before_intervention')
    def test_failed_recovery_finite_horizon(self):
        r=classify([0,1,2,3],[.1,0,0,0],actual_removal_day=1.5)
        self.assertEqual(r['classification'],'failed_recovery_by_observed_horizon')
    def test_initially_dry_is_not_lost_and_recovered(self):
        r=classify([0,1,2,3],[0,.1,.1,.1],actual_removal_day=.5)
        self.assertEqual(r['classification'],'initially_dry')
    def test_truncation(self):
        r=classify([0,1,2],[.1,0,0],actual_removal_day=1,expected_horizon_day=3)
        self.assertEqual(r['classification'],'unresolved_truncated')
    def test_interpolated_duration_bracket_uncertainty(self):
        r=classify([0,1,2,3],[.1,0,0,.1],actual_removal_day=1, persistence_days=.5)
        self.assertEqual(r['classification'],'unresolved_persistence_bracket')
    def test_requested_timestamp_never_substitutes_actual(self):
        r=classify([0,1,2,3,4],[.1,0,.1,.1,.1],actual_removal_day=2.5,requested_removal_day=1.5)
        self.assertEqual(r['removal_delay_days'],1)
        self.assertEqual(r['assessment_start_day'],2.5)
    def test_history_groups_ignore_future_interventions(self):
        base=dict(initial_state={'bed':[0,1]},physics={'h':200},boundary={'kind':'open'},
                  forcing_history=[{'time_day':0,'Q':.07}],event_day=2.48)
        a=history_key(**base,interventions=[])
        b=history_key(**base,interventions=[{'time_day':8,'clear':True}])
        c=history_key(**base,interventions=[{'time_day':12,'clear':True}])
        d=history_key(**base,interventions=[{'time_day':1,'clear':True}])
        self.assertEqual(a,b);self.assertEqual(b,c);self.assertNotEqual(a,d)
    def test_history_changes_physics(self):
        base=dict(initial_state={'bed':[0,1]},physics={'h':200},boundary={'kind':'open'},forcing_history=[{'time_day':0,'Q':.07}],event_day=2)
        a=history_key(**base);base['physics']['h']=100
        self.assertNotEqual(a,history_key(**base))
    def test_invalid(self):
        for t,f in [([0,0],[1,0]),([0,1],[1,float('nan')]),([0],[1])]:
            with self.assertRaises(ValueError):classify(t,f)
    def test_raw_adapter_requires_and_uses_actual_export(self):
        run=dict(removal_day=1.,days=4.,records=[{'day':t,'f_toe':f} for t,f in zip([0,1,2,3,4],[.1,0,.1,.1,.1])])
        with self.assertRaises(ValueError):classify_run(run)
        run['export']=[dict(t=1.5*86400)]
        result=classify_run(run)
        self.assertEqual(result['actual_removal_day'],1.5)
        self.assertEqual(result['removal_delay_days'],.5)


if __name__=='__main__':unittest.main()
