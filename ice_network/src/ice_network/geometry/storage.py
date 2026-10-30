"""Storage nodes: analytic stage-volume relation for a power-law basin.

Surface area at depth h above the basin bottom: A(h) = A_ref (h / h_ref)^q, q >= 0 (q = 0: prismatic).
Volume V(h) = A_ref h^(q+1) / ((q+1) h_ref^q); stage H = z_bottom + h. Both directions are closed-form, so no
interpolation table (and no interpolation error) enters the conservation audit.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def volume_from_depth(h, a_ref, h_ref, q):
    if h <= 0.0:
        return 0.0
    return a_ref * h ** (q + 1.0) / ((q + 1.0) * h_ref ** q)


@njit(cache=True)
def depth_from_volume(v, a_ref, h_ref, q):
    if v <= 0.0:
        return 0.0
    return ((q + 1.0) * h_ref ** q * v / a_ref) ** (1.0 / (q + 1.0))


@njit(cache=True)
def area_from_depth(h, a_ref, h_ref, q):
    return a_ref * (max(h, 0.0) / h_ref) ** q
