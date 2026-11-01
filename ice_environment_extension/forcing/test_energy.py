import unittest
from energy import (ColumnParameters, ColumnState, step_column, atmospheric_flux,
                    BareIceState, step_bare_ice)


class EnergyTests(unittest.TestCase):
    def test_bare_ice_exact_signed_sequence(self):
        a = step_bare_ice(BareIceState(100), 100, 10, 20, 5)
        self.assertEqual(a.state, BareIceState(0, 1.5))
        self.assertEqual(a.debris_storage_j_m2, 0)
        self.assertEqual(a.residual_j_m2, 0)
        b = step_bare_ice(a.state, 100, 10, -30, 5)
        self.assertEqual(b.state, BareIceState(250, 1.5))
        self.assertEqual(b.cold_storage_j_m2, -250)
        self.assertEqual(b.residual_j_m2, 0)

    def setUp(self):
        self.p = ColumnParameters(1000, 10, 334000)

    def test_cold_repaid_before_melt(self):
        a = step_column(ColumnState(0, 100), self.p, 10, 0, 5)
        self.assertEqual(a.state.cold_reserve_j_m2, 50)
        self.assertEqual(a.latent_j_m2, 0)
        b = step_column(a.state, self.p, 20, 0, 5)
        self.assertEqual(b.state.cold_reserve_j_m2, 0)
        self.assertEqual(b.latent_j_m2, 50)
        self.assertAlmostEqual(b.residual_j_m2, 0)

    def test_signed_cooling_and_long_step_stability(self):
        a = step_column(ColumnState(10, 0), self.p, 1e8, -100, 0)
        self.assertTrue(-10 <= a.state.debris_temperature_c <= 10)
        self.assertGreater(a.state.cold_reserve_j_m2, 0)
        self.assertAlmostEqual(a.residual_j_m2, 0, delta=1e-5)

    def test_storage_conservation_many_steps(self):
        s = ColumnState(3, 10000)
        total = 0
        for n in range(100):
            r = step_column(s, self.p, 13, 100 if n % 2 else -80, 3)
            self.assertAlmostEqual(r.residual_j_m2, 0, delta=1e-9)
            total += r.atmospheric_j_m2 + r.water_j_m2
            s = r.state
        stored = 1000 * (s.debris_temperature_c - 3) + 10000 - s.cold_reserve_j_m2 + 334000 * s.melted_kg_m2
        self.assertAlmostEqual(total, stored, delta=1e-8)

    def test_invalid_and_freezing(self):
        with self.assertRaises(ValueError):
            step_column(ColumnState(0, 0), self.p, 1, 0, -1)
        with self.assertRaises(NotImplementedError):
            step_column(ColumnState(0, 0), self.p, 1, 0, 0, refreeze=True)

    def test_atmosphere_signs(self):
        args = dict(face_shortwave_w_m2=100, incoming_longwave_w_m2=300,
                    surface_temperature_c=0, air_temperature_c=5,
                    air_specific_humidity_kg_kg=0, wind_speed_m_s=2,
                    pressure_pa=80000, albedo=.2, emissivity=.95,
                    sensible_exchange_coefficient=.002, moisture_exchange_coefficient=.002,
                    air_heat_capacity_j_kg_k=1005, latent_vaporization_j_kg=2.5e6,
                    dry_air_gas_constant_j_kg_k=287, stefan_boltzmann_w_m2_k4=5.67e-8)
        r = atmospheric_flux(**args)
        self.assertEqual(r['shortwave'], 80)
        self.assertGreater(r['sensible'], 0)
        self.assertLess(r['latent'], 0)
        self.assertAlmostEqual(r['net'], sum(r[k] for k in ('shortwave', 'longwave', 'sensible', 'latent')))
        with self.assertRaises(ValueError):
            atmospheric_flux(**{**args, 'surface_temperature_c': -1})


if __name__ == '__main__':
    unittest.main()
