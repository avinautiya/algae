import importlib.util
from pathlib import Path
import unittest

from .energy import Parameters, State, AreaContent, remap_area, step, step_bare


class EnergyTests(unittest.TestCase):
    def setUp(self):
        self.p = Parameters(200., 4., 334000.)

    def test_reference_constant_area(self):
        path = Path(__file__).resolve().parents[1] / 'vendor/reference_energy.py'
        spec = importlib.util.spec_from_file_location('reference_energy', path)
        ref = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ref)
        a = 7.
        old = ref.ColumnState(3., 80., 2.)
        expected = ref.step_column(old, ref.ColumnParameters(200., 4., 334000.), 20., 12., 3.)
        result = step(State(a, a*200*3, a*80, a*2*334000), self.p, 20., 12., 3.)
        self.assertAlmostEqual(result.state.debris_energy_j / (a*200), expected.state.debris_temperature_c)
        self.assertAlmostEqual(result.state.cold_debt_j/a, expected.state.cold_reserve_j_m2)
        self.assertAlmostEqual(result.state.cumulative_latent_j/(a*334000), expected.state.melted_kg_m2)
        self.assertAlmostEqual(result.residual_j, 0.)

    def test_cold_first_and_signed_ports(self):
        s = State(2., 0., 100.)
        zero = step(s, self.p, 10., 0., 0.)
        self.assertEqual(zero.state, s)
        warm = step(s, self.p, 10., 0., 3.)
        self.assertEqual(warm.state.cold_debt_j, 40.)
        self.assertEqual(warm.latent_j, 0.)
        cooled = step(warm.state, self.p, 10., 0., -3.)
        self.assertEqual(cooled.state.cold_debt_j, 100.)
        melt = step(s, self.p, 10., 0., 8.)
        self.assertEqual(melt.latent_j, 60.)

    def test_expansion_requires_declared_content(self):
        s = State(2., 100., 40., 20.)
        with self.assertRaises(ValueError):
            remap_area(s, 4.)
        r = remap_area(s, 4., incoming=AreaContent(0., 5., 0.))
        self.assertEqual(r.state.debris_energy_j, 100.)
        self.assertEqual(r.state.cold_debt_j, 50.)
        self.assertEqual(r.residual_j, 0.)
        self.assertEqual(r.state.accounted_energy_j-s.accounted_energy_j, -10.)

    def test_contraction_exports_all_content(self):
        s = State(4., 100., 40., 20.)
        r = remap_area(s, 2., outgoing=AreaContent(25., 10., 5.))
        self.assertEqual(r.state, State(2., 50., 20., 10.))
        self.assertEqual((r.debris_transfer_j, r.cold_debt_transfer_j, r.latent_transfer_j), (-50., -20., -10.))
        self.assertEqual(r.residual_j, 0.)
        with self.assertRaises(ValueError):
            remap_area(s, 2., outgoing=AreaContent(0., 30., 0.))

    def test_invalid_numeric(self):
        for bad in (float('nan'), float('inf'), True):
            with self.assertRaises(ValueError):
                State(1., bad, 0.)
            with self.assertRaises(ValueError):
                step(State(1., 0., 0.), self.p, 1., bad, 0.)
        with self.assertRaises(ValueError):
            step(State(1., 0., 0.), self.p, 1e308, 1e308, 0.)
        with self.assertRaises(ValueError):
            State(0., 0., 0.)
        with self.assertRaises(ValueError):
            Parameters(0., 4., 1.)

    def test_bare_ports(self):
        r = step_bare(State(2., 0., 100.), 10., 2., 4.)
        self.assertEqual(r.state, State(2., 0., 0., 20.))
        self.assertEqual(r.residual_j, 0.)
        r = step_bare(r.state, 10., -2., -4.)
        self.assertEqual(r.state, State(2., 0., 120., 20.))
        with self.assertRaises(ValueError):
            step_bare(State(1., 1., 0.), 1., 0., 0.)


if __name__ == '__main__':
    unittest.main()
