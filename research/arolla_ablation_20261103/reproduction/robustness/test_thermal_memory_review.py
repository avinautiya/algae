import unittest
from thermal_memory_independent_review import steady_single_cell,eigen_check,continuum_eigen_convergence,demo
class IndependentThermalReview(unittest.TestCase):
 def test_single_cell_two_half_cells(self):
  r=steady_single_cell();self.assertLess(r['max_flux_error_W_m2'],1e-10)
 def test_exact_discrete_eigenmode(self):
  r=eigen_check();self.assertLess(r['max_temperature_error_against_exact_discrete_eigenmode_C'],1e-10)
 def test_continuum_eigen_refinement(self):
  errors=[r['max_temperature_error_C'] for r in continuum_eigen_convergence()]
  self.assertLess(errors[1],errors[0]/3);self.assertLess(errors[2],errors[1]/3)
 def test_periodic_transfer_and_endpoint_storage(self):
  for d in [.18,.6]:
   r=demo(d)
   self.assertLess(r['relative_transfer_error'],.02)
   self.assertLess(abs(r['energy_relative_residual']),1e-10)
   mismatch=r['total_incoming_heat_J_m2']-r['total_basal_heat_J_m2']-r['final_minus_initial_storage_J_m2']
   self.assertLess(abs(mismatch),1e-4)
   self.assertLess(r['numerical_basal_harmonic_amplitude_W_m2'],r['quasistatic_amplitude_W_m2'])
if __name__=='__main__':unittest.main()
