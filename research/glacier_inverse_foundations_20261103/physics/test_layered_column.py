import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from paths import require_piramide
sys.path.insert(0,str(require_piramide()/'model'))
from cold_ice_column import simulate
from layered_column import predict


class LayeredPhysicsTests(unittest.TestCase):
    def test_homogeneous_matches_frozen_phase_model(self):
        Ts=[5.,-5.,10.];N=12;depth=.1;initial=[5*(1-(i+.5)/N) for i in range(N)]
        old=simulate(debris_conductivity_W_mK=1.,debris_heat_capacity_J_m3K=2e6,debris_depth_m=depth,
            initial_debris_temperature_C=initial,initial_ice_temperature_C=[0.]*24,
            surface_temperature_C=Ts,interval_seconds=900.,dt_seconds=30.,ice_depth_m=.2)
        new=predict(k_upper=1.,k_lower=1.,capacity=2e6,depth=depth,layer_depth=.037,
            initial_depths=[0.,depth],initial_temperatures=[5.,0.],surface_samples=Ts+[10.],sample_depths=[.025,.05,.075],
            interval_s=900.,dt_s=30.,debris_cells=N,ice_cells=24,ice_depth=.2)
        self.assertAlmostEqual(new['melt_mmwe'],old['melt_mmwe'],places=11)
        self.assertLess(abs(new['relative_energy_residual']),1e-12)

    def test_unaligned_two_layer_exact_series_steady_solution(self):
        depth=.1;split=.037;k1=2.;k2=.4;Ts=10.;q=Ts/(split/k1+(depth-split)/k2)
        def T(z):return Ts-q*(min(z,split)/k1+max(z-split,0)/k2)
        # Include actual discontinuity in the initial interpolation.
        r=predict(k_upper=k1,k_lower=k2,capacity=2e6,depth=depth,layer_depth=split,
            initial_depths=[0.,split,depth],initial_temperatures=[Ts,T(split),0.],surface_samples=[Ts,Ts],
            sample_depths=[.02,.06,.09],interval_s=900.,dt_s=30.,debris_cells=12,ice_cells=24,ice_depth=.2)
        self.assertAlmostEqual(r['melt_mmwe'],q*900/333700,places=11)
        self.assertLess(abs(r['relative_energy_residual']),1e-12)

    def test_zero_gradient_uniform_cold_limit_and_invalid_inputs(self):
        r=predict(k_upper=1.,k_lower=.5,capacity=2e6,depth=.1,layer_depth=.05,
            initial_depths=[0.,.1],initial_temperatures=[0.,0.],surface_samples=[0.,0.],sample_depths=[.05])
        self.assertEqual(r['melt_mmwe'],0.)
        self.assertEqual(r['predicted_temperature_C'],[[0.],[0.]])

    def test_two_layer_warm_temperature_invariance_does_not_identify_melt(self):
        def run(scale):
            return predict(k_upper=2.*scale,k_lower=.5*scale,capacity=2e6*scale,depth=.1,layer_depth=.037,
                initial_depths=[0.,.1],initial_temperatures=[0.,0.],surface_samples=[0.,4.,8.,3.],
                sample_depths=[.025,.05,.075],interval_s=900.,dt_s=30.,debris_cells=16,ice_cells=24,ice_depth=.2)
        a,b=run(1.),run(1.7)
        self.assertEqual(a['cold_interface_steps'],0)
        self.assertEqual(b['cold_interface_steps'],0)
        self.assertLess(max(abs(u-v) for x,y in zip(a['predicted_temperature_C'],b['predicted_temperature_C']) for u,v in zip(x,y)),1e-13)
        self.assertAlmostEqual(b['melt_mmwe']/a['melt_mmwe'],1.7,places=12)


if __name__=='__main__':unittest.main()
