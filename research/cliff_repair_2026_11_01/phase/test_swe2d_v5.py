import unittest
import numpy as np
from swe2d_v5 import Sim,RHO_I,RHO_W,CP,LF

def flat_cell_domain():
    s=Sim(dx=.25,Lpool=.5,Lr=.5,Ltail=.25,Wl=.25,Q=0.,Tw=.5,
          heff=200.,closed_out=True,ramp_ratio=0)
    s.wall[:]=False;s.zb[:]=0;s.ice[:]=False;s.ice[1:-1,1:-1]=True
    s.reset_state(.2)
    return s

class IntegratedSourceTests(unittest.TestCase):
    def test_closed_actual_kernel_mass_and_latent_heat(self):
        s=flat_cell_domain();v0=s.liquid_volume();e0=s.sensible_heat();h0=s.h.copy();z0=s.zb.copy()
        result=s.step_physical(2.)
        dv=s.liquid_volume()-v0;vi=result['source']['ice_volume_m3']
        self.assertGreater(vi,0.)
        self.assertAlmostEqual(RHO_W*dv,RHO_I*vi,delta=2e-12)
        self.assertAlmostEqual(e0-s.sensible_heat(),RHO_I*LF*vi,delta=2e-7)
        dz=z0-s.zb
        np.testing.assert_allclose((s.zb+s.h)[1:-1,1:-1],(z0+h0-(1-RHO_I/RHO_W)*dz)[1:-1,1:-1],atol=2e-16)
        self.assertLess(abs(result['ledgers']['physical_split_thermal_latent_residual_J']),2e-7)
    def test_source_preserves_momentum_sensible_thermal_dilution(self):
        s=flat_cell_domain();s._window(.1,True)
        h=s.h.copy();ht=s.hT.copy();s.hu[1:-1,1:-1]=.01;s.hv[1:-1,1:-1]=.02;mx=s.hu.copy();my=s.hv.copy()
        s.lower_bed(.1)
        np.testing.assert_array_equal(s.hT,ht);np.testing.assert_array_equal(s.hu,mx);np.testing.assert_array_equal(s.hv,my)
        self.assertTrue(np.all(s.h[1:-1,1:-1]>h[1:-1,1:-1]))
        self.assertTrue(np.all(s.hT[1:-1,1:-1]/s.h[1:-1,1:-1]<ht[1:-1,1:-1]/h[1:-1,1:-1]))
    def test_ghost_melt_is_excluded(self):
        s=flat_cell_domain();s._window(.1,True);s.ice[0,:]=True;s.last_qacc[0,:]=1e20
        before=s.zb[0].copy();h=s.h[0].copy();r=s.lower_bed(.1)
        np.testing.assert_array_equal(s.zb[0],before);np.testing.assert_array_equal(s.h[0],h)
        self.assertLess(r[0],1e-4)
    def test_no_repeated_source_without_fresh_heat(self):
        s=flat_cell_domain();s.step_physical(.1);before=s.h.copy()
        with self.assertRaises(RuntimeError):s.lower_bed(.1)
        np.testing.assert_array_equal(s.h,before)
    def test_explicit_reset_accounts_mass_and_heat(self):
        s=flat_cell_domain();s.step_physical(.1);s.reset_state(.3);r=s.combined_ledgers()
        self.assertLess(abs(r['executed_state_water_volume_residual_m3']),2e-15)
        self.assertLess(abs(r['executed_state_sensible_residual_J']),2e-7)
        self.assertGreater(abs(r['explicit_reset_water_m3']),0.)
        self.assertFalse(r['physical_split_only'])
        self.assertIsNone(r['physical_split_thermal_latent_residual_J'])
        with self.assertRaises(RuntimeError):s.lower_bed(.1)
    def test_quasisteady_heat_time_is_not_physical_energy_certificate(self):
        s=flat_cell_domain();s._window(.1,True);s.lower_bed(100.)
        self.assertAlmostEqual(s.last_phase_source['heat_time_scaling'],1000.)
        self.assertIsNone(s.combined_ledgers()['physical_split_thermal_latent_residual_J'])
        self.assertFalse(s.combined_ledgers()['physical_split_only'])
        self.assertGreater(abs(s.last_phase_source['latent_minus_window_heat_J']),0.)
    def test_zero_heat_limit(self):
        s=flat_cell_domain();before=s.zb.copy();r=s.step_physical(.2,cool=False)
        self.assertEqual(r['ice_melt_m3'],0.);np.testing.assert_array_equal(s.zb,before)
    def test_multistep_closed_cell_ledgers(self):
        s=flat_cell_domain()
        for _ in range(8):s.step_physical(.5)
        r=s.combined_ledgers();self.assertLess(abs(r['combined_ice_water_mass_residual_kg']),1e-10)
        self.assertLess(abs(r['physical_split_thermal_latent_residual_J']),2e-6)
        self.assertAlmostEqual(s.t,4.)
    def test_bad_interval_cannot_mutate(self):
        s=flat_cell_domain();before=s.h.copy()
        for dt in [0.,-1.,float('nan'),float('inf'),True]:
            with self.assertRaises(ValueError):s.step_physical(dt)
        np.testing.assert_array_equal(s.h,before)
    def test_uniform_counterfactual_preserves_total_and_region_api(self):
        s=flat_cell_domain();s.incision='uniform';mask=np.zeros_like(s.ice);mask[1:-1,1:-1]=True;s.regions={'all':mask}
        s._window(.1,True);p=s.region_heat()['all'];v,_,_=s.lower_bed(.1)
        self.assertAlmostEqual(v,p*.1/(RHO_I*LF),delta=1e-18)
        self.assertAlmostEqual(s.region_vol['all'],v);self.assertEqual(len(s.last_reach),2)
    def test_hydraulics_matches_existing_api_and_flags_stale_diagnostics(self):
        s=flat_cell_domain();r=s.hydraulics(tmax=.2,window=.1,min_t=.2)
        self.assertIn('converged',r);self.assertTrue(s.diagnostics()['hydraulic_window_state_is_current'])
        s.lower_bed(.1);self.assertFalse(s.diagnostics()['hydraulic_window_state_is_current'])
    def test_phase_birth_zero_initial_water_is_supported_by_source(self):
        s=flat_cell_domain();s.h[:]=0;s.hu[:]=0;s.hv[:]=0;s.hT[:]=0
        # Direct source fixture intentionally assumes an externally supplied heat
        # window; this is not a dry-cell water-only heat-transfer claim.
        s.last_qacc=np.zeros_like(s.zb);s.last_qacc[1,1]=100.
        s._heat_ready=True;s._heat_epoch=1;s._source_window=dict(elapsed_s=1.,ice_heat_removed_J=100.*s.dx**2)
        s.lower_bed(1.);self.assertGreater(s.h[1,1],0.);self.assertEqual(s.hT[1,1],0.)
    def test_open_actual_kernel_boundary_export_and_combined_ledgers(self):
        s=Sim(dx=.25,Lpool=.5,Lr=.5,Ltail=.25,Wl=.25,Q=.01,Tw=.5,ramp_ratio=0)
        for _ in range(10):s.step_physical(.2)
        r=s.combined_ledgers()
        self.assertLess(abs(r['combined_ice_water_mass_residual_kg']),2e-10)
        self.assertLess(abs(r['physical_split_thermal_latent_residual_J']),2e-6)
        self.assertGreater(abs(r['hydraulic_boundary_net_water_m3']),0.)
        self.assertGreater(r['generated_water_m3'],0.)
        self.assertIn('phase_to_absolute_section_throughflow',s.last_phase_source['route_source_diagnostics']['toe'])
    def test_settling_cannot_be_relabelled_as_physical_time(self):
        s=flat_cell_domain();s.hydraulics(tmax=.2,window=.1,min_t=.2);s.step_physical(.1)
        self.assertFalse(s.combined_ledgers()['physical_split_only'])
        self.assertIsNone(s.combined_ledgers()['physical_split_thermal_latent_residual_J'])
    def test_unsupported_constructor_inputs_fail_before_run(self):
        for cfg in [dict(Q=-1.),dict(Tw=-.1),dict(dx=0.),dict(heff=0.),dict(cfl=.9)]:
            with self.assertRaises(ValueError):Sim(**cfg)
    def test_covered_physical_sources_refused_before_heat_or_flow_mutation(self):
        for mode in ['solid','porous']:
            s=flat_cell_domain();height=np.zeros((s.nx,s.ny));height[1,0]=.01;s.set_deposit(height,mode=mode)
            before=[s.h.copy(),s.hT.copy(),s.zb.copy()]
            with self.assertRaises(ValueError):s.step_physical(.1)
            for old,new in zip(before,[s.h,s.hT,s.zb]):np.testing.assert_array_equal(old,new)
            self.assertEqual(len(s.hydraulic_window_log),0)
    def test_quasisteady_covered_source_assumption_is_explicit(self):
        s=flat_cell_domain();height=np.zeros((s.nx,s.ny));height[1,0]=.01;s.set_deposit(height)
        s._window(.1,True);s.lower_bed(.1)
        self.assertGreater(s.last_phase_source['covered_ice_melt_m3'],0.)
        self.assertIn('not physically verified',s.last_phase_source['under_debris_source_transport'])
    def test_geometry_mutation_requires_new_heat_window(self):
        s=flat_cell_domain();s._window(.1,True);s.remove_debris()
        with self.assertRaises(RuntimeError):s.lower_bed(.1)
    def test_uniform_counterfactual_not_native_physical_claim(self):
        s=flat_cell_domain();s.incision='uniform'
        with self.assertRaises(ValueError):s.step_physical(.1)

if __name__=='__main__':unittest.main()
