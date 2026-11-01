import math
import unittest
from standard_closures import WaterIceProperties, nusselt, integrate_constant_cells

P = WaterIceProperties(1.8e-6, 13.5, .56, 917., 334000.)

class Comparators(unittest.TestCase):
    def test_cooling_heating_ratio_exact(self):
        self.assertAlmostEqual(nusselt(1e5, 13.5, 'db_published_04') /
                               nusselt(1e5, 13.5, 'db_cooling_03'), 13.5**.1)

    def test_gnielinski_independent_rearrangement(self):
        re, pr = 120000., 13.5
        f = (1.82*math.log10(re)-1.64)**-2
        expected = (re-1000)*pr/(8/f+12.7*math.sqrt(8/f)*(pr**(2/3)-1))
        self.assertAlmostEqual(nusselt(re, pr, 'gnielinski_smooth'), expected)

    def test_integral_exact_overlap_and_zero_temperature(self):
        cells = [dict(start_s=a,end_s=b,temperature_C=t,velocity_m_s=.8,hydraulic_diameter_m=1.3)
                 for a,b,t in [(0,10,0),(10,20,.2)]]
        got = integrate_constant_cells(cells, 5, 15, P, 'db_published_04')
        re = .8*1.3/P.nu_m2_s
        expected = .023*re**.8*P.prandtl**.4*P.water_k_W_m_K/1.3*.2*5
        self.assertAlmostEqual(got['thermal_energy_J_m2'], expected)
        self.assertAlmostEqual(got['predicted_lowering_m']*P.ice_density_kg_m3*P.latent_J_kg, expected)

    def test_gaps_extrapolation_domains_and_negative_temperature_refused(self):
        with self.assertRaises(ValueError): nusselt(100, 13.5, 'db_published_04')
        with self.assertRaises(ValueError): nusselt(1e5, 13.5, 'gnielinski_supplied_f')
        cell=dict(start_s=0,end_s=10,temperature_C=.1,velocity_m_s=.8,hydraulic_diameter_m=1.3)
        with self.assertRaises(ValueError): integrate_constant_cells([cell],0,11,P,'db_published_04')
        with self.assertRaises(ValueError): integrate_constant_cells([dict(cell,temperature_C=-.1)],0,10,P,'db_published_04')

if __name__ == '__main__': unittest.main()
