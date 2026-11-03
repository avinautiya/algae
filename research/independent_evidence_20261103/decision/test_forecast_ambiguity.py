import unittest
from forecast_ambiguity import power,equivalent_temperature,demonstration
class Tests(unittest.TestCase):
 def test_training_equivalence_and_forecast_difference(self):
  d=demonstration()
  for p in d['pairs']:self.assertAlmostEqual(p['training_power_W'],100.,12)
  self.assertGreater(d['future_power_difference_W'],20.)
 def test_film_limit(self):
  self.assertAlmostEqual(power(1e12,1.,100.,1.)/100.,1.,9)
 def test_finite_supply_limit(self):self.assertLessEqual(power(10.,1.,10000.,.5),5.)
 def test_independent_temperature_resolves_h(self):
  F=1000.;A=1.;h=300.;T=.5;P=power(F,A,h,T)
  inferred=-F/A*__import__('math').log1p(-P/(F*T))
  self.assertAlmostEqual(inferred,h,12)
if __name__=='__main__':unittest.main()
