import itertools,random,unittest
from bounded_mae import bounded_loss_difference,issued_forecast_bound
class BoundTests(unittest.TestCase):
    def test_shared_reading_bound(self):
        r=issued_forecast_bound([2,3,4],[4,5,6],[2,3,4],.1)
        self.assertAlmostEqual(r['certified_conservative_difference_interval'][1],-1.6)
        self.assertEqual(r['conclusion'],'candidate_strictly_better_throughout_box')
    def test_corner_only_can_miss_internal_maximum(self):
        def loss(e):return (abs(1-(2+e))+abs(3-(2+e))-2*abs(2-(2+e)))/2
        self.assertEqual(loss(-2),0);self.assertEqual(loss(2),0);self.assertEqual(loss(0),1)
    def test_random_dependent_error_boxes(self):
        random.seed(20261103)
        for trial in range(100):
            a=[random.uniform(-3,3) for _ in range(3)];b=[random.uniform(-3,3)for _ in range(3)]
            A=[[random.uniform(-2,2)for _ in range(4)]for _ in range(3)];B=[[random.uniform(-2,2)for _ in range(4)]for _ in range(3)]
            ds=[random.random()for _ in range(4)];r=bounded_loss_difference(a,b,A,B,ds);lo,hi=r['certified_conservative_difference_interval']
            for _ in range(20):
                e=[random.uniform(-d,d)for d in ds]
                v=sum(abs(x+sum(c*z for c,z in zip(row,e)))-abs(y+sum(c*z for c,z in zip(other,e)))for x,y,row,other in zip(a,b,A,B))/3
                self.assertLessEqual(lo-1e-12,v);self.assertLessEqual(v,hi+1e-12)
if __name__=='__main__':unittest.main()
