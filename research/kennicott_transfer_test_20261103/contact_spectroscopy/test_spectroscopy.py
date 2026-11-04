import unittest,math,cmath
from .spectroscopy import transfer,OMEGA,harmonic_extract,bounded_fit,jacobian_rank

def fv_harmonic(alpha,b,omega,N,L=.17,x=.1):
 # Independent cell-centered finite-volume frequency solve, C=1,k=alpha.
 dx=L/N;g=alpha/dx;gd=alpha/(b+dx/2);diag=[2*g+1j*omega*dx]*N;diag[0]=3*g+1j*omega*dx;diag[-1]=g+gd+1j*omega*dx;rhs=[0j]*N;rhs[0]=2*g;lo=[-g]*(N-1);hi=[-g]*(N-1)
 for i in range(1,N):
  f=lo[i-1]/diag[i-1];diag[i]-=f*hi[i-1];rhs[i]-=f*rhs[i-1]
 sol=[0j]*N;sol[-1]=rhs[-1]/diag[-1]
 for i in range(N-2,-1,-1):sol[i]=(rhs[i]-hi[i]*sol[i+1])/diag[i]
 u=x/dx-.5;i=int(u);w=u-i;return sol[i]*(1-w)+sol[i+1]*w

def observations(alpha,b):
 return [{'omega_rad_s':w,'ratio_real':transfer(alpha,b,w).real,'ratio_imag':transfer(alpha,b,w).imag,'ratio_error_scale':.01} for w in [0,OMEGA,2*OMEGA]]
class Tests(unittest.TestCase):
 def test_DC_and_null_limits(self):
  self.assertAlmostEqual(transfer(1e-6,0,0).real,.07/.17);self.assertAlmostEqual(transfer(1e-6,.1,0).real,.17/.27)
  for w in [OMEGA,2*OMEGA]:self.assertAlmostEqual(abs(transfer(1e-6,0,w)-transfer(1e-6,1e-12,w)),0,places=9)
 def test_fv_independent_refinement(self):
  exact=transfer(8e-7,.035,OMEGA);errors=[]
  for n in [34,68,136]:errors.append(abs(fv_harmonic(8e-7,.035,OMEGA,n)-exact))
  self.assertGreater(errors[0]/errors[1],3);self.assertGreater(errors[1]/errors[2],3);self.assertLess(errors[-1],1e-5)
 def test_synthetic_recovery_and_null_comparator(self):
  obs=observations(8e-7,.035);fit=bounded_fit(obs);null=bounded_fit(obs,null=True)
  self.assertLess(abs(fit['alpha_m2_s']/8e-7-1),.002);self.assertLess(abs(fit['b_kR_m']-.035),.0002);self.assertLess(fit['weighted_objective'],null['weighted_objective'])
 def test_fit_independent_FV_ratios(self):
  obs=[]
  for w in [0,OMEGA,2*OMEGA]:
   z=fv_harmonic(1.1e-6,.06,w,256);obs.append({'omega_rad_s':w,'ratio_real':z.real,'ratio_imag':z.imag,'ratio_error_scale':.01})
  f=bounded_fit(obs);self.assertLess(abs(f['alpha_m2_s']/1.1e-6-1),.003);self.assertLess(abs(f['b_kR_m']-.06),.0003)
 def test_OLS_ratio_phase_and_trend(self):
  alpha,b=8e-7,.035;times=[i*600 for i in range(2017)];U0=10;u=[];y=[]
  for t in times:
   z1=complex(2,-1);z2=complex(1.2,.3);trend=.1*(t-604800)/86400
   u.append(U0+trend+(z1*cmath.exp(1j*OMEGA*t)).real+(z2*cmath.exp(2j*OMEGA*t)).real)
   y.append(U0*transfer(alpha,b,0).real+.03*(t-604800)/86400+(transfer(alpha,b,OMEGA)*z1*cmath.exp(1j*OMEGA*t)).real+(transfer(alpha,b,2*OMEGA)*z2*cmath.exp(2j*OMEGA*t)).real)
  f=harmonic_extract(times,u,y)
  for row,w in zip(f['observations'],[0,OMEGA,2*OMEGA]):self.assertLess(abs(complex(row['ratio_real'],row['ratio_imag'])-transfer(alpha,b,w)),1e-12)
 def test_rank_is_two_for_rich_signal(self):
  s=jacobian_rank(8e-7,.035,observations(8e-7,.035))['singular_values'];self.assertGreater(s[1],1)
 def test_b0_truth_not_contact_falsepositive(self):
  f=bounded_fit(observations(8e-7,0));self.assertLess(f['b_kR_m'],.0002)
 def test_bad_inputs(self):
  for a,b in [(-1,0),(1,-1),(float('nan'),0)]:
   with self.assertRaises(ValueError):transfer(a,b,OMEGA)
 def test_unresolved_input_refused(self):
  with self.assertRaises(ValueError):harmonic_extract(list(range(20)),[0.]*20,[0.]*20)
if __name__=='__main__':unittest.main()
