import unittest
from .spectroscopy import transfer,OMEGA
from .geometry_equivalence import shifted_boundary_transfer
class Tests(unittest.TestCase):
 def test_exact_complex_shift(self):
  for a,b,L,x in [(8e-7,.02,.17,.1),(7e-7,.013,.23,.225),(2e-7,.2,.4,.2)]:
   for w in [0,OMEGA,2*OMEGA]:self.assertLess(abs(transfer(a,b,w,L,x)-shifted_boundary_transfer(a,b,w,L,x)),2e-15)
 def test_small_contact_real_shift_is_cubic(self):
  a,L,x=8e-7,.17,.1
  def e(b):return abs(transfer(a,b,OMEGA,L,x)-transfer(a,0,OMEGA,L+b,x))
  self.assertGreater(e(.002)/e(.001),7.7);self.assertLess(e(.002)/e(.001),8.3)
if __name__=='__main__':unittest.main()
