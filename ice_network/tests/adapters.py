"""TEST-ONLY adapters that let the production storage/ODE/ledger machinery run the Codex limiting cases.

The linear and Torricelli edge laws are prescribed test constitutive laws, NOT Manning/GVF channels. They exercise node storage, the
ODE integration, junction accounting and the cumulative ledgers, not the production hydraulics. The production hydraulics law is never
modified, and no expected value is computed from the production solver.
"""
import numpy as np

from ice_network.evolution.network import Network


class PrescribedLawNetwork(Network):
    """Edges carry `law` = ("linear", G) -> Q = G * max(H_a - datum, 0) or ("sqrt", C) -> Q = C * sqrt(max(H_a - datum, 0)); flow a -> b only."""

    def edge_state(self, e, h_a, h_b, t_a, t_b=0.0):
        kind, coef = e["law"]
        head = max(h_a - e["datum"], 0.0)
        q = coef * head if kind == "linear" else coef * np.sqrt(head)
        n = e["n"]
        out = dict(Q=q, dir=1, h=np.zeros(n), hf=np.zeros(n - 1), choked=np.zeros(n, np.int64), Qabs=q, phi_in=0.0, phi_exit=0.0,
                   phi_seg=np.zeros(n - 1), phi_tot=0.0, q_ice=np.zeros(n - 1), wall_frac=np.zeros(n - 1), t_out=t_a, q_melt=0.0)
        return out


def make_reservoir(area, edges_laws, datum=0.0):
    lakes = [dict(z_bottom=datum, a_ref=area, h_ref=1.0, q=0.0)]
    outlets = [dict(h_fixed=datum - 1.0) for _ in edges_laws]
    edges = [dict(a=0, b=1 + i, z=np.array([datum, datum]), w=np.array([1.0, 1.0]), zs=np.array([datum + 1, datum + 1]), length=1.0,
                  n_man=0.04, ke=0.5, law=law, datum=datum) for i, law in enumerate(edges_laws)]
    return PrescribedLawNetwork(lakes, outlets, edges, melting=False)
