"""Conditional sensor-error bounds, not a new inverse theorem or field model."""
import math


def residual_radius(mass_row, deltas, adjacent_edges, duration, discrepancy=0.0):
    """Edges are (G, node_i, node_j), each G>=0; fixed geometry and C.

    Assumes pointwise sensor bounds and a nonnegative independently supplied
    integrated model discrepancy. The returned radius encloses any coherent
    sensor trajectory satisfying those bounds, including shared endpoints.
    """
    if len(mass_row) != len(deltas):
        raise ValueError('mass and uncertainty shapes differ')
    vals = list(mass_row)+list(deltas)+[duration, discrepancy]
    if not all(math.isfinite(x) for x in vals):
        raise ValueError('finite inputs required')
    if min(deltas, default=0) < 0 or duration < 0 or discrepancy < 0:
        raise ValueError('nonnegative uncertainty and time required')
    out = 2*sum(abs(m)*d for m,d in zip(mass_row,deltas))+discrepancy
    for conductance, i, j in adjacent_edges:
        if not math.isfinite(conductance) or conductance < 0:
            raise ValueError('passive finite conductance required')
        out += duration*conductance*(deltas[i]+deltas[j])
    return out


def integrated_residual(mass_row, initial, final, integrals, adjacent_edges):
    """Interior integrated energy residual in J/m²."""
    if not (len(mass_row)==len(initial)==len(final)==len(integrals)):
        raise ValueError('shape mismatch')
    return sum(m*(v-u) for m,u,v in zip(mass_row,initial,final))+sum(
        g*(integrals[i]-integrals[j]) for g,i,j in adjacent_edges)


def basal_energy_interval(conductance_low, conductance_high, observed_integral,
                          duration, delta):
    """Fixed continuously-melting basal temperature zero; no cold-ice claim.

    Returns signed heat bounds. A negative lower bound is NOT positive melt.
    Supports potentially negative integrals without silently clipping them.
    """
    vals=[conductance_low,conductance_high,observed_integral,duration,delta]
    if not all(math.isfinite(x) for x in vals):
        raise ValueError('finite inputs required')
    if conductance_low<0 or conductance_high<conductance_low or min(duration,delta)<0:
        raise ValueError('invalid bounds')
    lo=observed_integral-duration*delta
    hi=observed_integral+duration*delta
    corners=[g*x for g in (conductance_low,conductance_high) for x in (lo,hi)]
    return min(corners),max(corners)
