"""Synthetic algebra checks only; no camera data."""
import unittest
import numpy as np

def issued_ratio(air,theta,train,beta):
    reference=np.mean(theta[train]);dose=air*np.exp(beta*(theta-reference))
    return np.sum(dose[~train])/np.sum(dose[train])
class Checks(unittest.TestCase):
 def setUp(self):
    self.air=np.array([1.,2.,3.,2.,2.,3.]);self.theta=np.array([.07,.08,.06,.09,.10,.08]);self.train=np.array([True,True,True,False,False,False])
 def test_constant_additive_calibration_invariance(self):
    x=issued_ratio(self.air,self.theta,self.train,8.)
    for bias in [-.031,.02,.1]:self.assertAlmostEqual(x,issued_ratio(self.air,self.theta+bias,self.train,8.),13)
 def test_fixed_beta_drift_bound(self):
    beta=8.;delta=.031;x=issued_ratio(self.air,self.theta,self.train,beta)
    e=np.where(self.train,-delta,delta);y=issued_ratio(self.air,self.theta+e,self.train,beta)
    self.assertAlmostEqual(y/x,np.exp(2*beta*delta),13)
 def test_multiplicative_calibration_does_not_cancel(self):
    x=issued_ratio(self.air,self.theta,self.train,8.);y=issued_ratio(self.air,self.theta*1.3,self.train,8.)
    self.assertNotAlmostEqual(x,y,6)
if __name__=='__main__':unittest.main()
