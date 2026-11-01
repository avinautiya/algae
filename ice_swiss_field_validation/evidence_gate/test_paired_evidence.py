import unittest
from paired_evidence import compare, worst_case_advantage


class EvidenceTests(unittest.TestCase):
    def test_measurement_uncertainty_can_destroy_point_advantage(self):
        self.assertGreater(worst_case_advantage([0, 0], [0, 0], [1, 1]), 0)
        self.assertLess(worst_case_advantage([-2, 2], [0, 0], [1, 1]), 0)

    def test_prediction_uncertainty_is_not_ignored(self):
        self.assertEqual(worst_case_advantage([0, 0], [-2, 2], [1, 1]), -1)

    def test_extremal_formula_matches_dense_numerical_enumeration(self):
        for obs, cand, base in [([-3, 4], [-1, 2], [-2, 3]), ([0, 4], [1, 2], [-3, -1]), ([1, 2], [4, 5], [1.3, 1.7])]:
            brute = min(max(base[0]-y, 0, y-base[1])-max(abs(cand[0]-y), abs(cand[1]-y)) for y in [obs[0]+(obs[1]-obs[0])*i/10000 for i in range(10001)])
            self.assertAlmostEqual(worst_case_advantage(obs, cand, base), brute, places=3)

    def test_pooled_gain_does_not_establish_consistency(self):
        rows = [dict(id='a', history='A', target_type='bed_lowering', unit='m', observation=[0,0], predictions={'new':[0,0], 'old':[10,10]}),
                dict(id='b', history='B', target_type='bed_lowering', unit='m', observation=[0,0], predictions={'new':[1,1], 'old':[0,0]})]
        self.assertFalse(compare(rows, 'new', ['old'])['consistent_within_supplied_bounds'])

    def test_duplicate_target_is_rejected(self):
        row=dict(id='same', history='A', target_type='bed_lowering', unit='m', observation=[0,0], predictions={'new':[0,0], 'old':[1,1]})
        with self.assertRaises(ValueError): compare([row, row], 'new', ['old'])


if __name__ == '__main__':
    unittest.main()
