import math,unittest
from temperature_model import predict,conditional_melt


class NativeSignedColumnTests(unittest.TestCase):
    def solve(self,upper,lower,initial,**extra):
        opts=dict(domain_depth_m=.5,initial_depths_m=[0,.25,.5],initial_temperature_C=initial,
                  upper_temperature_C=upper,lower_temperature_C=lower,interval_seconds=1800,
                  sample_depths_m=[.125,.25,.375],cells=20,dt_seconds=90)
        opts.update(extra);return predict(5e-7,**opts)

    def test_signed_linear_exact_no_clipping(self):
        r=self.solve([-3.,-3.,-3.],[1.,1.,1.],[-3.,-1.,1.])
        for row in r['predicted_temperature_C']:
            for actual,expected in zip(row,[-2.,-1.,0.]):self.assertAlmostEqual(actual,expected,places=12)
        self.assertTrue(r['subzero_debris_computed_without_clipping'])
        self.assertFalse(r['melt_prediction'])
        self.assertLess(abs(r['normalized_energy_residual_K_m']),1e-12)

    def test_native_halfhour_reconstruction_explicit(self):
        left=self.solve([0.,20.,0.],[0.,0.,0.],[0.,0.,0.],hold='left')
        right=self.solve([0.,20.,0.],[0.,0.,0.],[0.,0.,0.],hold='right')
        self.assertEqual(left['predicted_temperature_C'][1],[0.,0.,0.])
        self.assertGreater(right['predicted_temperature_C'][1][0],0.)
        self.assertEqual(len(left['interval_basal_normalized_heat_K_m']),2)

    def test_bad_inputs_and_no_extrapolation(self):
        for options in [dict(dt_seconds=True),dict(cells=1),dict(hold='unknown'),
                        dict(sample_depths_m=[.6]),dict(upper_temperature_C=[1.,float('nan')]),
                        dict(initial_temperature_C=[1,0,0]),dict(lower_temperature_C=[0])]:
            with self.assertRaises(ValueError):self.solve([0.,0.],[0.,0.],[0.,0.,0.],**options)

    def test_conditional_melt_gate_and_units(self):
        contract=dict(base_is_ice_interface=True,basal_temperature_zero_supported=True,
                      pore_phase_change_negligible_supported=True,bulk_capacity_independent_of_test_melt=True,
                      source='synthetic fixture; not a field declaration')
        r=conditional_melt([.1,.2],bulk_heat_capacity_J_m3K=2e6,latent_heat_J_kg=334000,ice_interface_contract=contract)
        self.assertAlmostEqual(r['conditional_melt_mmwe'],600000/334000)
        for key in list(contract)[:-1]:
            bad=dict(contract);bad[key]=False
            with self.assertRaises(ValueError):conditional_melt([.1],bulk_heat_capacity_J_m3K=2e6,latent_heat_J_kg=334000,ice_interface_contract=bad)
        with self.assertRaises(ValueError):conditional_melt([-.1],bulk_heat_capacity_J_m3K=2e6,latent_heat_J_kg=334000,ice_interface_contract=contract)


if __name__=='__main__':unittest.main()
