import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from .solver import Column, SurfaceBoundary, ensemble, simulate


def source():
    return dict(status='declared_synthetic', source='analytic fixture',
                location='test_solver.py', method='declared exact test problem',
                uses_evaluation_melt=False)


def column(n=12, profile=None, depth=.1, k=1., capacity=2e6):
    sources={name:source() for name in ('conductivity','heat_capacity','depth','initial_profile','basal_ice_state')}
    sources['heat_capacity']['quantity']='bulk_debris_volumetric_heat_capacity'
    sources['basal_ice_state'].update(state='temperate_melting_interface',temperature_C=0.)
    return Column(k, capacity, depth,
                  tuple(profile if profile is not None else [0.]*n),sources)


def boundary(values):
    s = source(); s['quantity'] = 'debris_surface_temperature'
    return SurfaceBoundary(tuple(values), s)


def run(c, b, horizon=3600., step=60., **kwargs):
    return simulate(c, b, horizon_s=horizon, max_horizon_s=horizon,
                    dt_s=step, synthetic=True, **kwargs)


class ThermalMemoryTests(unittest.TestCase):
    def test_exact_steady_linear_solution_and_units(self):
        for n in (1, 2, 17):
            profile = [10*(1-(j+.5)/n) for j in range(n)]
            r = run(column(n, profile), boundary([10.]))
            self.assertLess(max(abs(a-b) for a,b in zip(profile,r['final_temperature_C'])),1e-11)
            self.assertAlmostEqual(r['basal_heat_J_m2'],3600*1.*10/.1,places=7)
            self.assertAlmostEqual(r['melt_mm_water_equivalent'],360000/333700,places=12)

    def test_discrete_eigenmode_backward_euler(self):
        n=19; step=10.; horizon=1000.; alpha=1/2e6; dz=.1/n
        initial=[math.sin(math.pi*(j+.5)/n) for j in range(n)]
        rate=4*alpha/dz**2*math.sin(math.pi/(2*n))**2
        factor=(1+step*rate)**(-int(horizon/step))
        r=run(column(n,initial),boundary([0.]),horizon,step)
        self.assertLess(max(abs(actual-factor*x) for actual,x in zip(r['final_temperature_C'],initial)),2e-13)

    def test_continuum_eigenmode_spatial_temporal_refinement(self):
        errors=[]
        for n in (8,16,32):
            initial=[math.sin(math.pi*(j+.5)/n) for j in range(n)]
            r=run(column(n,initial),boundary([0.]),1000.,5*(8/n)**2)
            factor=math.exp(-1/2e6*(math.pi/.1)**2*1000)
            errors.append(max(abs(actual-factor*x) for actual,x in zip(r['final_temperature_C'],initial)))
        self.assertGreater(errors[0]/errors[1],3.8)
        self.assertGreater(errors[1]/errors[2],3.9)
        self.assertLess(errors[-1],.0005)

    def test_first_order_temporal_refinement_against_discrete_mode(self):
        n=16;initial=[math.sin(math.pi*(j+.5)/n) for j in range(n)]
        rate=4*(1/2e6)/(.1/n)**2*math.sin(math.pi/(2*n))**2
        exact=math.exp(-rate*1000);errors=[]
        for step in (40.,20.,10.):
            r=run(column(n,initial),boundary([0.]),1000.,step)
            errors.append(max(abs(v-exact*x) for v,x in zip(r['final_temperature_C'],initial)))
        self.assertGreater(errors[0]/errors[1],1.98)
        self.assertGreater(errors[1]/errors[2],1.99)

    def test_exact_discrete_ledger_with_cooling_and_initial_storage(self):
        r=run(column(13,[4.]*13),boundary([12.,0.,3.]),10000.,733.)
        self.assertLess(abs(r['discrete_energy_relative_residual']),2e-13)
        for row in r['history']:
            change=row['storage_change_J_m2']
            flux=row['step_s']*(row['q_surface_W_m2']-row['q_basal_W_m2'])
            self.assertLess(abs(change-flux),1e-7)
        self.assertTrue(any(x['q_surface_W_m2']<0 for x in r['history']))

    def test_hourly_boundaries_and_partial_final_step(self):
        r=run(column(),boundary([1.,6.,2.]),10000.,733.)
        endpoints={x['time_s'] for x in r['history']}
        self.assertTrue({3600.,7200.,10000.}.issubset(endpoints))
        for row in r['history']:
            self.assertEqual(row['surface_temperature_C'],[1.,6.,2.][int((row['time_s']-row['step_s'])//3600)])

    def test_no_source_or_target_fitted_input_refused(self):
        for key in ('conductivity','heat_capacity','depth','initial_profile','basal_ice_state'):
            for bad in (None,{},dict(source(),status='target_fitted'),dict(source(),uses_evaluation_melt=True)):
                c=column();c.sources[key]=bad
                with self.assertRaises(ValueError):run(c,boundary([1.]))
        with self.assertRaises(ValueError):
            simulate(column(),boundary([1.]),horizon_s=3600,max_horizon_s=3600,dt_s=60)

    def test_bulk_capacity_and_temperate_interface_source_gates(self):
        c=column();c.sources['heat_capacity']['quantity']='rock_volumetric_heat_capacity'
        with self.assertRaises(ValueError):run(c,boundary([1.]))
        for value in (-1.,None,True,float('nan')):
            c=column();c.sources['basal_ice_state']['temperature_C']=value
            with self.assertRaises(ValueError):run(c,boundary([1.]))
        c=column();c.sources['basal_ice_state']['state']='cold_ice'
        with self.assertRaises(ValueError):run(c,boundary([1.]))

    def test_independent_declared_sources_permitted_without_certification(self):
        c=column(); b=boundary([1.])
        for s in list(c.sources.values())+[b.source]:s['status']='independent_literature'
        r=simulate(c,b,horizon_s=3600,max_horizon_s=3600,dt_s=60)
        self.assertFalse(r['synthetic'])
        self.assertFalse(r['validated_against_field_melt'])

    def test_nonfinite_boolean_negative_refreezing_inputs_refused(self):
        for value in (float('nan'),float('inf'),-1.,True):
            with self.assertRaises(ValueError):run(column(profile=[value]),boundary([1.]))
            with self.assertRaises(ValueError):run(column(),boundary([value]))
            for param in ('k','depth','capacity'):
                with self.assertRaises(ValueError):run(column(**{param:value}),boundary([1.]))
        for param in ('k','depth','capacity'):
            with self.assertRaises(ValueError):run(column(**{param:0.}),boundary([1.]))

    def test_timing_coverage_dt_and_limit_refused(self):
        for step in (0.,-1.,True,float('nan')):
            with self.assertRaises(ValueError):run(column(),boundary([1.]),step=step)
        with self.assertRaises(ValueError):run(column(),boundary([1.]),horizon=3601.)
        with self.assertRaises(ValueError):
            simulate(column(),boundary([1.]),horizon_s=3600,max_horizon_s=3599,dt_s=60,synthetic=True)
        b=boundary([1.]);b.source['quantity']='air_temperature'
        with self.assertRaises(ValueError):run(column(),b)
        with self.assertRaises(ValueError):run(column(1,k=8e306),boundary([1.]),step=1e-300)

    def test_thermal_memory_initial_state_and_delayed_flux(self):
        cold=run(column(),boundary([10.,0.]),7200.,30.)
        warm=run(column(profile=[4.]*12),boundary([10.,0.]),7200.,30.)
        self.assertGreater(warm['melt_mm_water_equivalent'],cold['melt_mm_water_equivalent'])
        self.assertGreater(sum(x['step_s']*x['q_basal_W_m2'] for x in cold['history'] if x['time_s']>3600),0.)

    def test_ensemble_bounds_only_finite_scenarios(self):
        r=ensemble([column(),column(profile=[3.]*12)],boundary([5.]),
                   horizon_s=3600,max_horizon_s=3600,dt_s=30,synthetic=True)
        values=[x['melt_mm_water_equivalent'] for x in r['scenarios']]
        self.assertEqual(r['melt_mm_water_equivalent_range'],[min(values),max(values)])
        self.assertFalse(r['continuous_parameter_bound_established'])
        self.assertFalse(r['statistical_confidence_interval'])
        with self.assertRaises(ValueError):ensemble([],boundary([1.]))

    def test_cli_synthetic_output_and_overwrite_refusal(self):
        root=Path(__file__).resolve().parents[1]
        c=column();b=boundary([2.])
        payload=dict(scenarios=[c.__dict__],boundary=b.__dict__,
                     options=dict(horizon_s=3600,max_horizon_s=3600,dt_s=120,synthetic=True))
        with tempfile.TemporaryDirectory() as directory:
            inp=Path(directory)/'input.json';out=Path(directory)/'out.json'
            inp.write_text(json.dumps(payload))
            command=[sys.executable,'-m','thermal_memory',str(inp),str(out)]
            completed=subprocess.run(command,cwd=root,capture_output=True,text=True)
            self.assertEqual(completed.returncode,0,completed.stderr)
            self.assertTrue(json.loads(out.read_text())['scenarios'][0]['melt_mm_water_equivalent']>0)
            self.assertNotEqual(subprocess.run(command,cwd=root,capture_output=True).returncode,0)


if __name__=='__main__':unittest.main()
