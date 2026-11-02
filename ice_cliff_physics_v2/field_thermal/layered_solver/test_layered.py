import unittest
import numpy as np
from layered import PreparedLayered
class Tests(unittest.TestCase):
 def test_homogeneous_sine(self):
  t=np.arange(101)*100.;z=np.linspace(0,.2,10001);a=1e-6;errors=[]
  for n in [20,40,80]:
   p=PreparedLayered(t,np.zeros(len(t)),np.zeros(len(t)),[.05,.1,.15],z,np.sin(np.pi*z/.2),cells=n)
   exact=np.exp(-a*(np.pi/.2)**2*t[:,None])*np.sin(np.pi*np.array([.05,.1,.15])/.2)
   errors.append(np.max(abs(p.predict(a,a)[1:]-exact[1:])))
  self.assertGreater(errors[0]/errors[1],3.8);self.assertGreater(errors[1]/errors[2],3.8)
 def test_series_resistance(self):
  z=np.linspace(0,.2,1001);a=1e-6;b=4e-6
  resistance=np.minimum(z,.1)/a+np.maximum(z-.1,0)/b
  T=3*(1-resistance/(.1/a+.1/b))
  p=PreparedLayered([0,300,600],[3]*3,[0]*3,[.05,.1,.15],z,T)
  np.testing.assert_allclose(p.predict(a,b),np.tile(np.interp([.05,.1,.15],z,T),(3,1)),atol=1e-12)
 def test_positivity_and_conservative_operator(self):
  p=PreparedLayered(np.arange(101)*300.,np.ones(101)*4,np.zeros(101),[.1],[0,.2],[4,0])
  d,e,f=p.operator(1e-7,3e-6);A=np.diag(d)+np.diag(e,1)+np.diag(e,-1)
  T=np.linspace(.2,2,p.n);top=4.;bot=0
  rhs=A@T;rhs[0]+=f[0]*top/p.dx**2;rhs[-1]+=f[-1]*bot/p.dx**2
  np.testing.assert_allclose(p.dx*rhs.sum(),(f[0]*(top-T[0])+f[-1]*(bot-T[-1]))/p.dx,atol=1e-15)
  H=p.history(1e-7,3e-6);self.assertGreaterEqual(H.min(),-1e-12);self.assertLessEqual(H.max(),4+1e-12)
 def test_invalid(self):
  with self.assertRaises(ValueError):PreparedLayered([0,600],[0,0],[0,0],[.1],[0,.2],[0,0])
if __name__=='__main__':unittest.main()
