import copy
import json
import math
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from ice_lake_coupling.kernel import (
    Constants, DomainError, LakeForcing, LakeState, PortTransfer, StageAreaTable,
    advance_lake, critical_discharge, exchange_port, thermal_control,
)
from ice_lake_coupling.simulation import run_case, validate_case

ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT/"ice_lake_coupling/cases/synthetic_lake.json"
C = Constants()


class StorageTests(unittest.TestCase):
    def test_exact_linear_area_integral(self):
        s = StageAreaTable((10, 12, 14), (100, 140, 180))
        for d in (0, 0.3, 1, 2, 2.7, 4):
            self.assertAlmostEqual(s.volume(10+d), 100*d+10*d*d, places=10)
            self.assertAlmostEqual(s.stage(100*d+10*d*d), 10+d, places=12)

    def test_random_piecewise_round_trip(self):
        rng = random.Random(622)
        # Area can locally decrease while volume still strictly increases.
        s = StageAreaTable((-1, 0, 2, 5), (50, 90, 65, 160))
        for _ in range(200):
            h = rng.uniform(-1, 5)
            self.assertAlmostEqual(s.stage(s.volume(h)), h, places=12)

    def test_no_extrapolation(self):
        s = StageAreaTable((0, 1), (10, 20))
        for value in (-0.01, 15.01):
            with self.assertRaises(DomainError):
                s.stage(value)
        with self.assertRaises(DomainError):
            s.volume(1.01)

    def test_bad_geometry_rejected(self):
        for heights, areas in (((0,0),(1,1)), ((0,1),(0,1)), ((0,1),(1,)), ((0,float("nan")),(1,1))):
            with self.subTest(heights=heights), self.assertRaises(DomainError):
                StageAreaTable(heights, areas)


class ConservationTests(unittest.TestCase):
    def setUp(self):
        self.storage = StageAreaTable((0, 10), (100, 100))

    def test_closed_lake_heating_exact(self):
        state = LakeState.from_temperature(500, 0.2)
        step = advance_lake(self.storage, state, LakeForcing(external_power_W=300), 100, lambda h: 0.0)
        self.assertEqual(step.state.volume_m3, 500)
        self.assertAlmostEqual(step.state.sensible_energy_J, state.sensible_energy_J+30000)
        self.assertEqual(step.water_residual_m3, 0)
        self.assertEqual(step.energy_residual_J, 0)

    def test_empty_lake_fills_and_mixes(self):
        step = advance_lake(self.storage, LakeState(0,0), LakeForcing(0.2,0.5), 100, lambda h: 0.0)
        self.assertEqual(step.state.volume_m3, 20)
        self.assertAlmostEqual(step.state.temperature(), 0.5)

    def test_empty_dry_and_dry_heating(self):
        step = advance_lake(self.storage, LakeState(0,0), LakeForcing(), 10, lambda h: 0.0)
        self.assertEqual(step.state, LakeState(0,0))
        with self.assertRaises(DomainError):
            advance_lake(self.storage, LakeState(0,0), LakeForcing(external_power_W=1), 10, lambda h: 0.0)

    def test_freezing_and_capacity_fail_explicitly(self):
        with self.assertRaises(DomainError):
            advance_lake(self.storage, LakeState.from_temperature(1,0.1), LakeForcing(external_power_W=-1e8), 10, lambda h: 0.0)
        with self.assertRaises(DomainError):
            advance_lake(self.storage, LakeState(1000,0), LakeForcing(1,0), 1, lambda h: 0.0)

    def test_drain_does_not_change_well_mixed_temperature(self):
        state = LakeState.from_temperature(500,0.37)
        for _ in range(100):
            step = advance_lake(self.storage, state, LakeForcing(), 50, lambda h: 0.1*h)
            self.assertLess(abs(step.energy_residual_J), 1e-6)
            state = step.state
        self.assertAlmostEqual(state.temperature(), 0.37, places=11)

    def test_stiff_drain_positive_without_volume_floor(self):
        step = advance_lake(self.storage, LakeState.from_temperature(500,0.4), LakeForcing(), 1000, lambda h: 1000*h)
        exact_be = 500/(1+1000*10)
        self.assertGreater(step.state.volume_m3, 0)
        self.assertAlmostEqual(step.state.volume_m3, exact_be, places=8)
        self.assertAlmostEqual(step.state.temperature(), 0.4, places=9)

    def test_external_port_antisymmetry_and_backflow(self):
        lake = LakeState.from_temperature(100,0.4)
        channel_volume, channel_energy = 7.0, C.rc*7*0.2
        for port in (PortTransfer(3,C.rc*3*0.4), PortTransfer(-2,-C.rc*2*0.2)):
            new = exchange_port(lake,port)
            self.assertEqual(new.volume_m3 + channel_volume + port.volume_m3, lake.volume_m3+channel_volume)
            self.assertAlmostEqual(new.sensible_energy_J+channel_energy+port.sensible_energy_J,
                                   lake.sensible_energy_J+channel_energy, places=6)

    def test_impossible_port_rejected(self):
        for volume, heat in ((1,-1), (0,2), (-1,2)):
            with self.assertRaises(DomainError):
                PortTransfer(volume,heat)
        with self.assertRaises(DomainError):
            exchange_port(LakeState(1,0),PortTransfer(2,0))

    def test_thermal_control_mass_energy_and_geometry(self):
        port = PortTransfer(12,C.rc*12*0.5)
        heat = thermal_control(port,60,2.5,1500)
        self.assertAlmostEqual(heat.latent_energy_J+heat.outflow.sensible_energy_J,port.sensible_energy_J,places=6)
        self.assertAlmostEqual(C.rho_water*heat.meltwater_m3,C.rho_ice*heat.ice_loss_m3,places=12)
        self.assertAlmostEqual(heat.lowering_m*2.5,heat.ice_loss_m3,places=12)
        self.assertEqual(heat.outflow.volume_m3,port.volume_m3+heat.meltwater_m3)

    def test_zero_flow_and_high_heat_transfer_limit(self):
        self.assertEqual(thermal_control(PortTransfer(0,0),1,1,100).latent_energy_J,0)
        port = PortTransfer(1,C.rc)
        heat = thermal_control(port,1,1,1e20)
        self.assertEqual(heat.latent_energy_J,port.sensible_energy_J)
        self.assertEqual(heat.outflow.sensible_energy_J,0)
        self.assertEqual(thermal_control(port,1,1,1e20,False).outflow,port)


class AnalyticalTests(unittest.TestCase):
    def test_critical_control_relation(self):
        depth = 0.2
        width = 0.7
        q = critical_discharge(1.5*depth,0,width)
        self.assertAlmostEqual(q/(width*depth),math.sqrt(C.gravity*depth),places=12)
        self.assertEqual(critical_discharge(-1,0,width),0)

    def test_nonlinear_drain_converges_to_analytic_solution(self):
        area, sill, h0, duration, width = 100,1,0.4,1000,0.2
        storage = StageAreaTable((0,3),(area,area))
        k = width*math.sqrt(C.gravity)*(2/3)**1.5
        exact = sill + (h0**(-0.5)+k*duration/(2*area))**(-2)
        errors=[]
        for dt in (20,10,5):
            state = LakeState.from_temperature(area*(sill+h0),0.2)
            for _ in range(int(duration/dt)):
                state = advance_lake(storage,state,LakeForcing(),dt,
                    lambda h: critical_discharge(h,sill,width)).state
            errors.append(abs(storage.stage(state.volume_m3)-exact))
        self.assertGreater(errors[0]/errors[1],1.8)
        self.assertGreater(errors[1]/errors[2],1.8)
        self.assertLess(errors[-1],0.001)

    def test_inflow_replacement_temperature_converges(self):
        storage=StageAreaTable((0,3),(100,100))
        exact=0.8+(0.2-0.8)*math.exp(-0.1*1000/100)
        errors=[]
        for dt in (20,10,5):
            state=LakeState.from_temperature(100,0.2)
            for _ in range(int(1000/dt)):
                # Q=0.1 at H=1, with the required zero-storage limit.
                state=advance_lake(storage,state,LakeForcing(0.1,0.8),dt,lambda h:0.1*h).state
            errors.append(abs(state.temperature()-exact))
        self.assertGreater(errors[0]/errors[1],1.8)
        self.assertGreater(errors[1]/errors[2],1.8)


class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.case=json.loads(CASE.read_text())

    def test_whole_case_ledgers_and_forcing_event_alignment(self):
        rows, summary=run_case(self.case,True)
        self.assertLess(summary["max_water_residual_relative"],1e-9)
        self.assertLess(summary["max_energy_residual_relative"],1e-12)
        self.assertLess(abs(summary["geometry_latent_volume_residual_m3"]),1e-11)
        self.assertEqual(rows[-1]["time_s"],21600)
        self.assertTrue({7200,14400,21600}.issubset({r["time_s"] for r in rows}))
        self.assertGreater(summary["sill_lowering_m"],0)
        self.assertTrue(all(r["lake_temperature_K"]>=0 and r["lake_volume_m3"]>=0 for r in rows))

    def test_fixed_control_no_melting(self):
        _, summary=run_case(self.case,False)
        self.assertEqual(summary["sill_lowering_m"],0)
        self.assertEqual(summary["totals"]["latent_J"],0)
        self.assertEqual(summary["totals"]["ice_loss_m3"],0)

    def test_geometry_substeps_are_conservative(self):
        self.case["max_lowering_per_step_m"]=0.00002
        rows, summary=run_case(self.case,True)
        self.assertGreater(summary["geometry_step_rejections"],0)
        self.assertLess(summary["max_water_residual_relative"],1e-9)
        self.assertTrue(all(a["sill_m"]-b["sill_m"]<=0.000020001 for a,b in zip(rows,rows[1:])))

    def test_unknown_parameter_and_bad_interval_rejected(self):
        self.case["control"]["hydrofracture"]=True
        with self.assertRaises(DomainError):
            validate_case(self.case)
        del self.case["control"]["hydrofracture"]
        self.case["forcing"][-1]["until_s"]=21599
        with self.assertRaises(DomainError):
            validate_case(self.case)

    def test_real_cli_runs_and_protects_existing_results(self):
        case=copy.deepcopy(self.case)
        case["duration_s"]=60
        case["forcing"]=case["forcing"][:1]
        case["forcing"][0]["until_s"]=60
        with tempfile.TemporaryDirectory() as temp:
            inp=Path(temp)/"case.json"; inp.write_text(json.dumps(case))
            out=Path(temp)/"out"
            command=[sys.executable,"-m","ice_lake_coupling","--case",str(inp),"--output",str(out)]
            first=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
            self.assertEqual(first.returncode,0,first.stderr)
            result=(out/"comparison.json").read_bytes()
            second=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
            self.assertNotEqual(second.returncode,0)
            self.assertEqual(result,(out/"comparison.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
