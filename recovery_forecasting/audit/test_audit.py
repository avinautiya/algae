import copy
import unittest
from .evidence import readiness,GATES
from .evaluate import evaluate,event_score


def ready():
    return dict(study_id='synthetic',claim_scope='conditional_computational',
                freeze=dict(configuration_sha256='a'*64,commit='b'*40,predictions_frozen_before_test=True),
                gates={g:dict(status='pass',record_sha256='c'*64) for g in GATES},
                campaign=dict(expected_case_ids=['a','b'],completed_case_ids=['b','a'],unresolved_failure_ids=[]))


def evaluation():
    return dict(target=dict(toe_fraction_threshold=.001,persistence_s=86400,after_ordinary_flow_return=True),
                models=['geometry','debris'],cases=[
                    dict(case_id='a',recovered=True,event=dict(kind='interval',lower_s=100,upper_s=120),
                         predictions=dict(geometry=dict(recovered=True,event_time_s=110),debris=dict(recovered=False,event_time_s=None))),
                    dict(case_id='b',recovered=False,event=dict(kind='right_censored',lower_s=1000),
                         predictions=dict(geometry=dict(recovered=False,event_time_s=None),debris=dict(recovered=True,event_time_s=500)))])


class AuditTests(unittest.TestCase):
    def test_conditional_metadata_complete(self):self.assertTrue(readiness(ready())['checklist_complete'])
    def test_actual_partition_failure_blocks(self):
        r=ready();r['gates']['hydraulic_partition']['status']='fail'
        self.assertFalse(readiness(r)['checklist_complete'])
    def test_missing_every_gate_blocks(self):
        for g in GATES:
            r=ready();del r['gates'][g]
            self.assertFalse(readiness(r)['checklist_complete'])
    def test_environment_requires_external(self):
        r=ready();r['claim_scope']='environmental_prediction'
        self.assertFalse(readiness(r)['checklist_complete'])
        r['external_validation']=dict(independent=True,status='pass',record_sha256='d'*64)
        self.assertTrue(readiness(r)['checklist_complete'])
    def test_campaign_missing_or_duplicate(self):
        for ids in (['a'],['a','a']):
            r=ready();r['campaign']['completed_case_ids']=ids
            self.assertFalse(readiness(r)['checklist_complete'])
    def test_nonfinite_rejected(self):
        r=ready();r['unknown']=float('nan')
        with self.assertRaises(ValueError):readiness(r)
    def test_classification_and_pairing(self):
        out=evaluate(evaluation())
        self.assertEqual(out['scores']['geometry']['accuracy_on_nonabstained'],1.)
        self.assertEqual(out['paired_comparisons'][0]['a_better'],2)
        self.assertEqual(out['scores']['geometry']['mean_distance_to_observed_event_bracket_s'],0.)
        self.assertTrue(out['scores']['debris']['event_scores'][1]['early_prediction_violation'])
    def test_interval_distances_not_exact_error(self):
        obs=dict(kind='interval',lower_s=10,upper_s=20)
        self.assertEqual(event_score(obs,5)['interval_distance_s'],5)
        self.assertEqual(event_score(obs,15)['interval_distance_s'],0)
        self.assertEqual(event_score(obs,25)['interval_distance_s'],5)
        self.assertTrue(event_score(obs,None)['missed_observed_event'])
    def test_censor_does_not_invent_mae(self):
        obs=dict(kind='right_censored',lower_s=10)
        self.assertFalse(event_score(obs,15)['early_prediction_violation'])
        self.assertNotIn('interval_distance_s',event_score(obs,15))
        self.assertFalse(event_score(obs,None)['early_prediction_violation'])
    def test_abstention_and_unresolved_not_counted_as_wrong(self):
        r=evaluation();r['cases'][0]['predictions']['geometry']['recovered']=None;r['cases'][1]['recovered']=None
        out=evaluate(r)
        self.assertEqual(out['scores']['geometry']['scored_classifications'],0)
        self.assertEqual(out['scores']['geometry']['abstentions'],1)
        self.assertEqual(out['scores']['geometry']['unresolved_labels'],1)
        self.assertIsNone(out['paired_comparisons'][0]['accuracy_difference'])
    def test_duplicate_or_missing_models_rejected(self):
        r=evaluation();r['cases'].append(copy.deepcopy(r['cases'][0]))
        with self.assertRaises(ValueError):evaluate(r)
        r=evaluation();del r['cases'][0]['predictions']['debris']
        with self.assertRaises(ValueError):evaluate(r)
    def test_invalid_target_and_event(self):
        r=evaluation();r['target']['after_ordinary_flow_return']=False
        with self.assertRaises(ValueError):evaluate(r)
        with self.assertRaises(ValueError):event_score(dict(kind='interval',lower_s=20,upper_s=10),15)
        with self.assertRaises(ValueError):event_score(dict(kind='right_censored',lower_s=0),float('nan'))


if __name__=='__main__':unittest.main()
