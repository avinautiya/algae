import itertools,math,random,sys,unittest
from pathlib import Path
from audit_actual_batch import audit
sys.path.insert(0,str(Path(__file__).parents[1]/'analysis'))
from bounded_mae import issued_forecast_bound
class IndependentBatchEvidenceTests(unittest.TestCase):
 def test_release_refused_before_any_actual_artifact_read(self):
  with self.assertRaisesRegex(ValueError,'authorize future evaluation release'):audit('/no/artifacts/here')
 def test_shared_endpoint_box_bound_contains_dense_independent_values(self):
  rng=random.Random(78)
  for _ in range(10):
   p=[rng.uniform(0,30) for _ in range(3)];q=[rng.uniform(0,30) for _ in range(3)];obs=[rng.uniform(0,30) for _ in range(3)];delta=2.
   lo,hi=issued_forecast_bound(p,q,obs,delta)['certified_conservative_difference_interval']
   for e in itertools.product([-delta,-delta/2,0.,delta/2,delta],repeat=4):
    truth=[y+e[i+1]-e[i] for i,y in enumerate(obs)]
    difference=math.fsum(abs(a-y)-abs(b-y) for a,b,y in zip(p,q,truth))/3
    self.assertGreaterEqual(difference+1e-12,lo);self.assertLessEqual(difference-1e-12,hi)
if __name__=='__main__':unittest.main()
