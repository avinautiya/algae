import json,math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'model'))
from cold_ice_column import simulate as frozen_solve
from contact_column import simulate


def inputs(*,R=.1,Ts=None,Nd=12,Ni=24,Td=None,Ti=None,step=30.,cadence=900.):
    return dict(debris_conductivity_W_mK=1.,debris_heat_capacity_J_m3K=2e6,debris_depth_m=.1,
        initial_debris_temperature_C=[0.]*Nd if Td is None else Td,
        initial_ice_temperature_C=[0.]*Ni if Ti is None else Ti,
        surface_temperature_C=[10.] if Ts is None else Ts,interval_seconds=cadence,dt_seconds=step,
        ice_depth_m=.2,contact_resistance_m2K_W=R)


class ContactColumnTests(unittest.TestCase):
    def check_accounting(self,r):
        self.assertLess(abs(r['relative_energy_residual']),2e-11)
        self.assertEqual(r['complementarity_max_product_C_kg_m2'],0.)
        self.assertGreaterEqual(r['contact_entropy_generation_J_m2K'],0.)
        self.assertGreaterEqual(r['physical_entropy_generation_J_m2K'],0.)
        self.assertGreaterEqual(r['backward_euler_entropy_defect_J_m2K'],0.)
        scale=max(abs(r['sensible_entropy_change_J_m2K']),abs(r['latent_entropy_J_m2K']),
                  abs(r['reservoir_entropy_in_J_m2K']),r['physical_entropy_generation_J_m2K'],1.)
        self.assertLess(abs(r['entropy_balance_residual_J_m2K'])/scale,2e-11)
        for row in r['history']:
            self.assertGreaterEqual(row['contact_entropy_generation_W_m2K'],0.)
            self.assertGreaterEqual(row['backward_euler_entropy_defect_J_m2K'],0.)
            if row['contact_law_residual_C'] is not None:self.assertLess(abs(row['contact_law_residual_C']),1e-11)
            self.assertLessEqual(row['interface_temperature_C'],0.)
            self.assertGreaterEqual(row['melt_kg_m2'],0.)

    def test_R_zero_reproduces_frozen_model_warm_and_cold(self):
        for temps in ([10.,-6.,4.],[-5.,-5.],[3.,8.,12.]):
            p=inputs(R=0.,Ts=temps);actual=simulate(**p);p.pop('contact_resistance_m2K_W');old=frozen_solve(**p)
            for key in ('melt_mmwe','surface_heat_J_m2','deep_heat_out_J_m2','final_sensible_storage_J_m2'):
                self.assertAlmostEqual(actual[key],old[key],places=10)
            self.assertEqual(actual['final_debris_temperature_C'],old['final_debris_temperature_C'])
            self.assertEqual(actual['final_ice_temperature_C'],old['final_ice_temperature_C'])
            self.check_accounting(actual)

    def test_exact_steady_warm_contact_series_resistance(self):
        for n in (1,7,24):
            R=.3;Ts=10.;q=Ts/(.1+R);profile=[Ts-q*(j+.5)*.1/n for j in range(n)]
            r=simulate(**inputs(R=R,Ts=[Ts],Nd=n,Td=profile,step=73.))
            self.assertAlmostEqual(r['melt_mmwe'],q*900/333700,places=12)
            self.assertAlmostEqual(r['history'][-1]['debris_contact_face_temperature_C'],q*R,places=12)
            self.assertEqual(r['history'][-1]['interface_temperature_C'],0.)
            self.assertLess(max(abs(a-b) for a,b in zip(profile,r['final_debris_temperature_C'])),1e-12)
            self.check_accounting(r)

    def test_exact_cold_series_resistance_and_entropy(self):
        R=.2;Ts=-5.;q=Ts/(.1+R+.2/2.1);face=Ts-q*.1;Tb=face-q*R
        Td=[Ts-q*(j+.5)*.1/12 for j in range(12)]
        Ti=[Tb*(1-(j+.5)/24) for j in range(24)]
        r=simulate(**inputs(R=R,Ts=[Ts],Td=Td,Ti=Ti,step=73.))
        self.assertEqual(r['melt_mmwe'],0.)
        self.assertAlmostEqual(r['history'][-1]['interface_temperature_C'],Tb,places=12)
        self.assertAlmostEqual(r['surface_heat_J_m2'],q*900,places=7)
        expected=900*R*q*q/((face+273.15)*(Tb+273.15))
        self.assertAlmostEqual(r['contact_entropy_generation_J_m2K'],expected,places=11)
        self.check_accounting(r)

    def test_infinite_resistance_is_exactly_insulated_and_serializable(self):
        r=simulate(**inputs(R=float('inf'),Ts=[10.,-5.],Ti=[-1.]*24))
        self.assertEqual(r['melt_mmwe'],0.)
        self.assertTrue(all(row['debris_interface_heat_W_m2']==0 for row in r['history']))
        self.assertEqual(r['contact_entropy_generation_J_m2K'],0.)
        json.dumps(r,allow_nan=False)
        self.check_accounting(r)

    def test_large_resistance_approaches_insulated_limit(self):
        a=simulate(**inputs(R=0.));b=simulate(**inputs(R=1e8))
        self.assertLess(b['melt_mmwe'],a['melt_mmwe']*1e-7)
        self.check_accounting(b)

    def test_transient_heating_cooling_entropy_and_phase(self):
        r=simulate(**inputs(R=.03,Ts=[10.,-8.,-8.,12.],step=17.))
        self.assertTrue(any(row['branch']=='cold_interface' for row in r['history']))
        self.assertTrue(any(row['melt_kg_m2']>0 for row in r['history']))
        self.check_accounting(r)

    def test_invalid_resistance_absolute_temperature_refused(self):
        for value in (-1.,float('nan'),True):
            with self.assertRaises(ValueError):simulate(**inputs(R=value))
        for value in (-273.15,-300.):
            with self.assertRaises(ValueError):simulate(**inputs(Ts=[value]))


if __name__=='__main__':unittest.main()
