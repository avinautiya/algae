import numpy as np
from ice_network.evolution.network import Network
from ice_network.hydraulics.channel import edge_flow


def chan(a, b, z0, drop, w, L=100.0, n=9, zs_off=1.0, nm=0.04, ke=0.5, msub=16):
    z = z0 - drop * np.linspace(0, 1, n)
    return dict(a=a, b=b, z=z, w=np.full(n, w), zs=z + zs_off, length=L, n_man=nm, ke=ke, msub=msub)


def lake(bottom=99.0, area=1000.0):
    return dict(z_bottom=bottom, a_ref=area, h_ref=1.0, q=0.0)


def q_of_stage(e, h_a, h_b):
    n = e["z"].size
    h = np.zeros(n); hf = np.zeros(n - 1); ch = np.zeros(n, np.int64)
    return edge_flow(h_a, h_b, e["z"], e["w"], e["n_man"], e["length"] / (n - 1), e["ke"], e.get("msub", 16), h, hf, ch)[0]


def const(v):
    v = np.atleast_1d(np.asarray(v, float))
    return lambda t: v
