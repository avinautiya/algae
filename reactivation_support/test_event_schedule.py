import unittest
from .event_schedule import schedule

class ScheduleTests(unittest.TestCase):
    def test_nonmultiple_horizon_completed_exactly(self):
        s=schedule(100.,30.)
        self.assertEqual(s[-1][1],100.)
        self.assertEqual(sum(x[2] for x in s),100.)
    def test_removal_not_delayed_to_next_regular_update(self):
        e=2.8*86400;s=schedule(15*86400,6*3600,[e])
        self.assertIn(e,[v[0] for v in s]+[s[-1][1]])
        self.assertTrue(all(dt<=6*3600 for _,_,dt in s))
    def test_floating_tick_alias_does_not_create_tiny_step(self):
        s=schedule(1.,.1,[.3])
        self.assertEqual(len(s),10)
        self.assertIn(.3,[x[0] for x in s])
    def test_invalid_event_rejected(self):
        with self.assertRaises(ValueError):schedule(1.,.1,[2.])
    def test_boundary_events_do_not_duplicate_intervals(self):
        s=schedule(1.,.5,[0.,1.])
        self.assertEqual(len(s),2)
if __name__=='__main__':unittest.main()
