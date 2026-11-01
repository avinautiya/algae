import unittest
import itertools
import numpy as np
from .bounds import opening_bound, exclusion


def connects(mask):
    nr,nc=mask.shape
    seen={(0,j) for j in range(nc) if mask[0,j]};front=list(seen)
    while front:
        i,j=front.pop()
        if i==nr-1:return True
        for ii,jj in ((i-1,j),(i+1,j),(i,j-1),(i,j+1)):
            if 0<=ii<nr and 0<=jj<nc and mask[ii,jj] and (ii,jj) not in seen:
                seen.add((ii,jj));front.append((ii,jj))
    return False


class BoundsTests(unittest.TestCase):
    def bound(self,z,**kw):return opening_bound(z,level_upper_m=0,dx_m=1,dy_m=1,**kw)

    def test_path_exact_against_exhaustive_subsets(self):
        rng=np.random.default_rng(71)
        for _ in range(12):
            z=rng.uniform(-.3,1,(3,3));best=float('inf')
            for bits in itertools.product((False,True),repeat=9):
                mask=np.array(bits).reshape(3,3)
                if connects(mask):best=min(best,float(np.maximum(z,0)[mask].sum()))
            self.assertAlmostEqual(self.bound(z)['path_only_volume_m3'],best)

    def test_width_relaxation_bounds_all_connected_wide_subsets(self):
        z=np.array([[.1,.7,.9],[.8,.1,.8],[.9,.7,.1]])
        bound=self.bound(z,minimum_cross_section_width_m=2,resolved_cell_route=True)
        best=float('inf')
        for bits in itertools.product((False,True),repeat=9):
            mask=np.array(bits).reshape(3,3)
            if (mask.sum(axis=1)>=2).all() and connects(mask):
                best=min(best,float(z[mask].sum()))
        self.assertLessEqual(bound['minimum_ice_volume_lower_bound_m3'],best+1e-12)
        self.assertGreaterEqual(bound['minimum_ice_volume_lower_bound_m3'],bound['path_only_volume_m3'])

    def test_partial_width_and_scaling(self):
        b=self.bound(np.ones((3,3)),minimum_cross_section_width_m=1.5)
        self.assertAlmostEqual(b['relaxed_width_volume_m3'],4.5)
        self.assertAlmostEqual(b['latent_energy_lower_bound_J'],4.5*917*334000)

    def test_level_ordering(self):
        z=np.array([[1.,2.],[3.,1.]])
        values=[opening_bound(z,level_upper_m=h,dx_m=1,dy_m=1)['latent_energy_lower_bound_J'] for h in (0,1,2,3)]
        self.assertTrue(all(a>=b for a,b in zip(values,values[1:])))

    def test_wall_and_zero_cost(self):
        z=np.zeros((3,3));wall=np.zeros_like(z,dtype=bool);wall[1,:]=True
        self.assertTrue(np.isinf(self.bound(z,blocked=wall,resolved_cell_route=True)['latent_energy_lower_bound_J']))
        self.assertEqual(self.bound(z)['latent_energy_lower_bound_J'],0)

    def test_strict_exclusion(self):
        b=self.bound([[1.]],resolved_cell_route=True)
        self.assertTrue(exclusion(b,0)['reconnection_excluded'])
        self.assertFalse(exclusion(b,b['latent_energy_lower_bound_J'])['reconnection_excluded'])

    def test_subcell_corridor_does_not_inherit_full_cell_cost(self):
        b=self.bound(np.ones((3,3)),minimum_cross_section_width_m=.1)
        self.assertAlmostEqual(b['minimum_ice_volume_lower_bound_m3'],.3)
        self.assertAlmostEqual(b['path_only_volume_m3'],3.)
        self.assertEqual(self.bound(np.ones((3,3)))['minimum_ice_volume_lower_bound_m3'],0.)

    def test_invalid(self):
        for z in ([[np.nan]],[],[[np.inf]]):
            with self.assertRaises(ValueError):self.bound(z)
        with self.assertRaises(ValueError):self.bound([[0]],minimum_cross_section_width_m=-1)
        with self.assertRaises(ValueError):self.bound([[0]],rho_i=True)
        with self.assertRaises(ValueError):exclusion({'latent_energy_lower_bound_J':float('nan')},1)


if __name__=='__main__':unittest.main()
