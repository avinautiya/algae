import math,unittest
from cold_ice_column import simulate
from temperature_model import predict


def run(Ts,*,Td=None,Ti=None,Nd=10,Ni=40,k=1.,C=2e6,depth=.1,ice_depth=2.,step=30.,interval=900.,profiles=False):
    return simulate(debris_conductivity_W_mK=k,debris_heat_capacity_J_m3K=C,debris_depth_m=depth,
        initial_debris_temperature_C=Td if Td is not None else [0.]*Nd,
        initial_ice_temperature_C=Ti if Ti is not None else [0.]*Ni,
        surface_temperature_C=Ts,interval_seconds=interval,dt_seconds=step,
        ice_depth_m=ice_depth,record_profiles=profiles)


class ColdIceComplementarityTests(unittest.TestCase):
    def check_ledger(self,r):
        self.assertLess(abs(r['relative_energy_residual']),3e-11)
        self.assertEqual(r['complementarity_max_product_C_kg_m2'],0.)
        for x in r['history']:
            self.assertLessEqual(x['interface_temperature_C'],0.)
            self.assertGreaterEqual(x['melt_kg_m2'],0.)
            self.assertLess(abs(x['energy_residual_J_m2']),1e-6)
            self.assertAlmostEqual(x['debris_interface_heat_W_m2']-x['ice_interface_heat_W_m2'],
                                   x['latent_heat_J_m2']/x['step_s'],places=8)

    def test_steady_warm_linear_exact_and_units(self):
        n=13;profile=[10*(1-(i+.5)/n) for i in range(n)]
        r=run([10.,10.],Td=profile,Ni=11,step=73.)
        self.assertAlmostEqual(r['melt_mmwe'],100*1800/333700,places=12)
        self.assertLess(max(abs(a-b) for a,b in zip(profile,r['final_debris_temperature_C'])),1e-12)
        self.assertEqual(r['final_ice_temperature_C'],[0.]*11)
        self.check_ledger(r)

    def test_steady_cold_two_material_resistances(self):
        Nd=10;Ni=20;depth=.1;ice_depth=.2;Ts=-5.;q=Ts/(depth+ice_depth/2.1);Tb=Ts-q*depth
        Td=[Ts-q*(i+.5)*depth/Nd for i in range(Nd)]
        Ti=[Tb*(1-(i+.5)/Ni) for i in range(Ni)]
        r=run([Ts],Td=Td,Ti=Ti,depth=depth,ice_depth=ice_depth,step=73.)
        self.assertEqual(r['melt_mmwe'],0.)
        self.assertAlmostEqual(r['history'][-1]['interface_temperature_C'],Tb,places=12)
        self.assertAlmostEqual(r['surface_heat_J_m2'],q*900,places=7)
        self.assertAlmostEqual(r['deep_heat_out_J_m2'],q*900,places=7)
        self.check_ledger(r)

    def test_cooling_stores_cold_content_not_negative_melt_or_water(self):
        r=run([-5.]*4)
        self.assertEqual(r['melt_mmwe'],0.)
        self.assertLess(r['final_sensible_storage_J_m2'],0.)
        self.assertLess(min(r['final_ice_temperature_C']),0.)
        self.assertFalse(r['refreezing_water_created'])
        self.assertTrue(all(t<=0 for t in r['final_ice_temperature_C']))
        self.check_ledger(r)

    def test_cold_ice_rewarming_delays_melt(self):
        cold=run([10.]*16,Ti=[-5.]*80,Ni=80,depth=.05,step=15.)
        warm=run([10.]*16,Ti=[0.]*80,Ni=80,depth=.05,step=15.)
        self.assertTrue(any(x['branch']=='cold_interface' for x in cold['history']))
        self.assertTrue(any(x['melt_kg_m2']>0 for x in cold['history']))
        first_cold=next(x['time_s'] for x in cold['history'] if x['melt_kg_m2']>0)
        first_warm=next(x['time_s'] for x in warm['history'] if x['melt_kg_m2']>0)
        self.assertGreater(first_cold,first_warm)
        self.assertLess(cold['melt_mmwe'],warm['melt_mmwe'])
        self.check_ledger(cold)

    def test_positive_warm_case_matches_previous_normalized_solver(self):
        N=16;k=1.3;C=1.6e6;cadence=1800.;dt=60.;Ts=[5.,12.,3.,8.]
        old=predict(k/C,domain_depth_m=.1,initial_depths_m=[0.,.1],initial_temperature_C=[5.,0.],
            upper_temperature_C=Ts+[8.],lower_temperature_C=[0.]*(len(Ts)+1),interval_seconds=cadence,
            sample_depths_m=[.05],cells=N,dt_seconds=dt)
        profile=[5*(1-(i+.5)/N) for i in range(N)]
        new=run(Ts,Td=profile,Nd=N,k=k,C=C,interval=cadence,step=dt)
        expected=C*math.fsum(old['interval_basal_normalized_heat_K_m'])/333700
        self.assertAlmostEqual(new['melt_mmwe'],expected,places=11)
        self.check_ledger(new)

    def test_spatial_refinement_cold_interface(self):
        values=[]
        for n in [8,16,32]:
            r=run([-5.,-5.],Nd=n,Ni=2*n,ice_depth=.2,step=1.)
            values.append(r['history'][-1]['interface_temperature_C']);self.check_ledger(r)
        self.assertLess(abs(values[2]-values[1]),abs(values[1]-values[0])*.4)

    def test_time_refinement_cold_interface(self):
        values=[]
        for step in [60.,30.,15.]:
            r=run([-5.,-5.],Nd=32,Ni=64,ice_depth=.2,step=step)
            values.append(r['history'][-1]['interface_temperature_C']);self.check_ledger(r)
        self.assertLess(abs(values[2]-values[1]),abs(values[1]-values[0])*.6)

    def test_nonfinite_bool_positive_ice_and_missing_profiles_refused(self):
        for bad in [float('nan'),float('inf'),True]:
            with self.assertRaises(ValueError):run([bad])
        with self.assertRaises(ValueError):run([0.],Ti=[1.])
        with self.assertRaises(ValueError):run([0.],Td=[])
        with self.assertRaises(ValueError):run([0.],k=0.)
        with self.assertRaises(ValueError):run([0.],step=0.)


if __name__=='__main__':unittest.main()
