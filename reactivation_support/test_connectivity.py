import unittest
from collections import deque
import numpy as np
from .connectivity import minimax_crest


def flood_fill_connects(z,level):
    nr,nc=z.shape;seen={(0,j) for j in range(nc) if z[0,j]<=level};queue=deque(seen)
    while queue:
        i,j=queue.popleft()
        if i==nr-1:return True
        for a,b in ((i-1,j),(i+1,j),(i,j-1),(i,j+1)):
            if 0<=a<nr and 0<=b<nc and z[a,b]<=level and (a,b) not in seen:seen.add((a,b));queue.append((a,b))
    return False

class ConnectivityTests(unittest.TestCase):
    def test_winding_corridor_beats_every_straight_column(self):
        z=np.array([[0,1,1],[0,0,1],[1,0,0],[1,1,0]],dtype=float)
        r=minimax_crest(z)
        self.assertEqual(r['crest'],0);self.assertEqual(r['straight_column_crest'],1)
    def test_uniform_transverse_barrier(self):
        z=np.zeros((5,4));z[2,:]=.2
        self.assertEqual(minimax_crest(z)['crest'],.2)
    def test_negative_datum_preserved(self):
        self.assertEqual(minimax_crest(np.full((4,3),-.08))['crest'],-.08)
    def test_diagonal_corner_contact_is_not_open_face(self):
        z=np.zeros((3,3));walls=np.ones((3,3),bool);np.fill_diagonal(walls,False)
        self.assertIsNone(minimax_crest(z,walls)['crest'])
    def test_independent_threshold_flood_fill(self):
        rng=np.random.default_rng(20261101)
        for _ in range(50):
            z=rng.uniform(-.1,.3,(4,5));r=minimax_crest(z)
            expected=next(v for v in sorted(z.ravel()) if flood_fill_connects(z,v))
            self.assertEqual(r['crest'],expected)
            self.assertEqual(max(z[tuple(v)] for v in r['path']),expected)
    def test_path_adjacency(self):
        r=minimax_crest(np.zeros((5,4)))
        self.assertEqual(r['path'][0][0],0);self.assertEqual(r['path'][-1][0],4)
        for a,b in zip(r['path'][:-1],r['path'][1:]):self.assertEqual(abs(a[0]-b[0])+abs(a[1]-b[1]),1)
    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):minimax_crest([[0,float('nan')],[0,0]])
if __name__=='__main__':unittest.main()
