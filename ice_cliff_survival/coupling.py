"""Explicitly UNFINISHED couplings (kept as hard failures so no success-shaped output can be produced).
Implemented and verified in milestone 1: 2-D heat conduction + moving ice interface (cfd/stefan2DFoam) with STATIC debris shielding (a prescribed conductance per interface row, e.g. from geometry.shield_conductance).
NOT implemented: mobile-grain transport coupled to the flow/heat/interface state (export capacity feeding back on the shield, deposition changing the covered rows in time), stream flow, free-surface hydraulics, burial dynamics."""
from . import UnsupportedRun
class MobileGrainCoupling:
    """placeholder for the dynamic shield: grain stock -> covered rows -> shielded interface conductance -> melt -> grain release. Every method raises."""
    reason = "mobile-grain transport coupling is unfinished: the shield in the 2-D solver is static; there is no verified closure for clast export/deposition driven by the flow, and no flow-heat-grain coupled verification"
    def __init__(self, *a, **k): raise UnsupportedRun(self.reason)
    def step(self, *a, **k): raise UnsupportedRun(self.reason)
def static_shield_rows(ny, Ly_face, z_f, theta_deg, beta_deg, k_deb):
    """hShield list for stefan2DFoam from a STATIC debris wedge: row j sits at height z_j = (j+0.5)/ny * Ly_face above the toe measured along the face coordinate -> vertical height z = s*sin(theta); k_deb [W m-1 K-1] is an assumed input."""
    import math
    from .geometry import shield_conductance
    return [shield_conductance(((j + 0.5) / ny * Ly_face) * math.sin(math.radians(theta_deg)), z_f, theta_deg, beta_deg, k_deb) for j in range(ny)]
