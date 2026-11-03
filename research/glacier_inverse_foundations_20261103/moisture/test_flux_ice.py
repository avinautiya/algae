import math,unittest
from flux_ice import simulate,soil_storage_components,prepare_flux,retrieve_from_plane


def run(q,*,ice=None,A=0.,Gamma=0.,step=30.,interval=900.,depth=2.):
    return simulate(interface_flux_W_m2=q,initial_ice_temperature_C=[0.]*64 if ice is None else ice,
        interval_seconds=interval,dt_seconds=step,ice_depth_m=depth,
        soil_interface_areal_capacity_J_m2K=A,initial_interface_temperature_C=Gamma)


class FluxIceTests(unittest.TestCase):
    def check(self,r):
        self.assertLess(abs(r['relative_energy_residual']),2e-11)
        self.assertEqual(r['complementarity_max_product'],0.)
        self.assertFalse(r['signed_flux_clipped'])
        self.assertFalse(r['refreezing_water_created'])
        for row in r['history']:
            self.assertLessEqual(row['interface_temperature_C'],0.)
            self.assertGreaterEqual(row['melt_kg_m2'],0.)
            self.assertLess(abs(row['energy_residual_J_m2']),1e-6)

    def test_positive_steady_exact_latent_units_with_soil_capacity(self):
        for A in (0.,20000.):
            r=run([25.]*4,A=A)
            self.assertAlmostEqual(r['melt_mmwe'],25*3600/333700,places=12)
            self.assertEqual(r['final_ice_temperature_C'],[0.]*64)
            self.check(r)

    def test_negative_flux_stores_cold_content_without_negative_melt(self):
        for A in (0.,20000.):
            r=run([-25.]*4,A=A)
            self.assertEqual(r['melt_mmwe'],0.)
            self.assertLess(r['final_interface_temperature_C'],0.)
            self.assertLess(r['final_sensible_storage_J_m2'],0.)
            self.assertTrue(any(x<0 for x in r['final_ice_temperature_C']))
            self.check(r)

    def test_single_cell_neumann_independent_analytic_oracle(self):
        q=-10.;h=90.;k=2.1;depth=.2;S=1925700.*depth;G=2*k/depth
        r=run([q],ice=[0.],step=h,interval=h,depth=depth)
        expected_T=q/(S/h+G);expected_Gamma=expected_T+q/G
        self.assertAlmostEqual(r['final_ice_temperature_C'][0],expected_T,places=14)
        self.assertAlmostEqual(r['final_interface_temperature_C'],expected_Gamma,places=14)
        self.check(r)

    def test_sheet_ice_two_node_independent_inverse(self):
        q=-10.;h=90.;depth=.2;A=20000.;S=1925700*depth;G=2*2.1/depth
        aa=A/h+G;dd=S/h+2*G;det=aa*dd-G*G
        expected_Gamma=q*dd/det;expected_T=q*G/det
        r=run([q],ice=[0.],A=A,step=h,interval=h,depth=depth)
        self.assertAlmostEqual(r['final_interface_temperature_C'],expected_Gamma,places=14)
        self.assertAlmostEqual(r['final_ice_temperature_C'][0],expected_T,places=14)
        self.check(r)

    def test_cooling_then_reheating_has_a_sensible_delay(self):
        r=run([-30.]*4+[50.]*8,A=20000.)
        first=next(x['time_s'] for x in r['history'] if x['melt_kg_m2']>0)
        self.assertGreater(first,3600.)
        self.assertLess(r['melt_mmwe'],50*7200/333700)
        self.assertTrue(any(x['branch']=='cold_interface_capacity' for x in r['history']))
        self.check(r)

    def test_linear_soil_integral_and_unknown_interface_term_once(self):
        C=2e6;z4=.375;z5=.5;d=.52;T4=2.;T5=.4;Gamma=-.3
        r=soil_storage_components(T4,T5,z4_m=z4,z5_m=z5,ice_depth_m=d,soil_capacity_J_m3K=C)
        p=(z4+z5)/2;Tp=(T4+T5)/2
        quadrature=C*((z5-p)*(Tp+T5)/2+(d-z5)*(T5+Gamma)/2)
        self.assertAlmostEqual(r['observed_storage_F_J_m2']+r['interface_areal_capacity_J_m2K']*Gamma,quadrature,places=8)
        self.assertAlmostEqual(r['interface_areal_capacity_J_m2K'],C*(d-z5)/2)

    def test_storage_subtraction_and_full_plane_ledger(self):
        F=[100000.,130000.,110000.]
        forcing=prepare_flux([30.,-10.],F,interval_seconds=1800.)
        self.assertAlmostEqual(forcing['effective_interface_flux_W_m2'][0],30.-30000/1800)
        r=retrieve_from_plane(plane_flux_W_m2=[30.,-10.],observed_storage_F_J_m2=F,
            initial_ice_temperature_C=[0.]*32,interval_seconds=1800.,dt_seconds=30.,soil_interface_areal_capacity_J_m2K=20000.)
        self.assertLess(abs(r['full_plane_energy_residual_J_m2']),1e-6)
        self.assertEqual(r['observed_soil_storage_change_J_m2'],10000.)
        self.check(r)

    def test_neumann_semiinfinite_cooling_spatial_refinement(self):
        q=-10.;t=1800.;k=2.1;alpha=k/1925700.;exact=2*q/k*math.sqrt(alpha*t/math.pi)
        errors=[]
        for n in (64,128,256):
            r=run([q],ice=[0.]*n,step=1.,interval=t)
            errors.append(abs(r['final_interface_temperature_C']-exact));self.check(r)
        self.assertLess(errors[1],errors[0]*.5)
        self.assertLess(errors[2],errors[1]*.6)

    def test_missing_nonfinite_and_phase_invalid_inputs_refused(self):
        for q in ([float('nan')],[True],[]):
            with self.assertRaises(ValueError):run(q)
        for options in (dict(A=-1.),dict(ice=[1.]),dict(Gamma=1.),dict(step=0.)):
            with self.assertRaises(ValueError):run([10.],**options)
        with self.assertRaises(ValueError):prepare_flux([1.,2.],[0.],interval_seconds=1800.)
        with self.assertRaises(ValueError):soil_storage_components(1.,0.,z4_m=.5,z5_m=.4,ice_depth_m=.52,soil_capacity_J_m3K=2e6)


if __name__=='__main__':unittest.main()
