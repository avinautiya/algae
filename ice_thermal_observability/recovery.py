"""Sharp finite-volume prediction ranges with independently evaluated LP duals.

These are deterministic bounds over a declared profile/measurement class, not
confidence intervals or certified enclosures of the continuum PDE.
"""
import numpy as np
from scipy.optimize import linprog


def extremum(objective, observation, data, uncertainty, lower, upper):
    """Maximize a scalar functional; return an attaining profile and dual bound."""
    c = np.asarray(objective, dtype=float)
    e = np.asarray(observation, dtype=float)
    d, eps = np.asarray(data, dtype=float), np.asarray(uncertainty, dtype=float)
    lo, hi = np.asarray(lower, dtype=float), np.asarray(upper, dtype=float)
    n = len(c)
    if e.ndim != 2 or e.shape[1] != n or e.shape[0] != len(d) or eps.shape != d.shape:
        raise ValueError("Incompatible functional/observation dimensions")
    if lo.shape != c.shape or hi.shape != c.shape:
        raise ValueError("One lower/upper bound per initial-profile cell required")
    if not all(np.isfinite(v).all() for v in (c,e,d,eps,lo,hi)):
        raise ValueError("Finite arrays required")
    if np.any(lo < 0) or np.any(lo > hi) or np.any(eps < 0):
        raise ValueError("Invalid temperature/measurement bounds")
    a = np.vstack([e,-e])
    b = np.concatenate([d+eps,-d+eps])
    result = linprog(-c, A_ub=a if len(d) else None, b_ub=b if len(d) else None,
                     bounds=list(zip(lo,hi)), method="highs",
                     options={"primal_feasibility_tolerance":1e-9,
                              "dual_feasibility_tolerance":1e-9})
    if not result.success:
        raise ValueError("Thermal information class infeasible/unresolved: "+result.message)
    profile = result.x
    if len(d):
        multipliers = -np.asarray(result.ineqlin.marginals)
        lam = multipliers[:len(d)]-multipliers[len(d):]
    else:
        lam = np.zeros(0)
    residual = c-e.T@lam
    dual = float(lam@d+np.abs(lam)@eps+np.maximum(residual,0)@hi+np.minimum(residual,0)@lo)
    primal = float(c@profile)
    feasibility = float(max(0,np.max(np.abs(e@profile-d)-eps))) if len(d) else 0.0
    gap = dual-primal
    scale = max(np.sum(np.abs(c)*np.maximum(np.abs(lo),np.abs(hi))),1e-12)
    if feasibility > 2e-8 or abs(gap) > 2e-8*scale:
        raise ArithmeticError("LP certificate did not meet validation tolerance")
    return dict(value=primal, variables=profile.tolist(), dual_upper_bound=dual,
                dual_gap=gap, dual_observation_weights=lam.tolist(),
                max_measurement_violation_K=feasibility,
                dual_residual_kernel=residual.tolist())


def sharp_interval(objective, observation, data, uncertainty, lower, upper):
    high = extremum(objective, observation, data, uncertainty, lower, upper)
    negative_low = extremum(-np.asarray(objective), observation, data, uncertainty, lower, upper)
    low = -negative_low["value"]
    upper_bound = high["value"]
    raw_low,raw_high=low,upper_bound
    scale=max(float(np.sum(np.abs(objective)*np.maximum(np.abs(lower),np.abs(upper)))),1e-12)
    resolution=2e-8*scale
    reordered=False
    if upper_bound<low:
        if low-upper_bound>resolution:
            raise ArithmeticError("Extremal values have inconsistent ordering")
        # Near-collapse LPs can disagree at roundoff. Preserve the raw result and
        # flag it, rather than returning a negative worst-case error radius.
        low,upper_bound=upper_bound,low
        reordered=True
    return dict(lower=low, upper=upper_bound, minimax=(low+upper_bound)/2,
                worst_case_absolute_error=(upper_bound-low)/2,
                raw_extremal_values=[raw_low,raw_high],extrema_reordered_for_roundoff=reordered,
                numerically_unresolved_width=upper_bound-low<=resolution,
                numerical_comparison_scale=resolution,
                lower_witness_variables=negative_low["variables"], upper_witness_variables=high["variables"],
                upper_certificate=high, negative_lower_certificate=negative_low)


def identifiable(objective, observation):
    """Unconstrained linear identifiability; boundaries may collapse a class too."""
    c, e = np.asarray(objective), np.asarray(observation)
    if not len(e):
        return dict(relative_rowspace_residual=1.0, rank=0)
    weights = np.linalg.lstsq(e.T,c,rcond=1e-12)[0]
    return dict(relative_rowspace_residual=float(np.linalg.norm(c-e.T@weights)/max(np.linalg.norm(c),1e-300)),
                rank=int(np.linalg.matrix_rank(e,tol=1e-12*np.linalg.norm(e))),
                reconstruction_weights=weights.tolist())
