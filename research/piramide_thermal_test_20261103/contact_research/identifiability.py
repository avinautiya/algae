"""Analytic model invariants and conditional necessary bounds, no field fitting."""
import math


def invariants(k,C,R,depth):
    values=[k,C,R,depth]
    if any(type(x) not in (int,float) or not math.isfinite(x) for x in values) or k<=0 or C<=0 or depth<=0 or R<0:
        raise ValueError('positive finite k,C,depth and nonnegative finite R required')
    return dict(diffusivity_m2_s=k/C,robin_length_m=k*R,
                steady_effective_depth_m=depth+k*R,
                temperature_only_absolute_heat_scale_identified=False)


def steady_temperate_sensor_bound(*,sensor_temperature_C,sensor_absolute_bound_C,mean_heat_flux_W_m2,
                                  temperate_steady_assumption_explicit):
    """Necessary constant-R bound for a sensor at/above the debris-side contact.

    q>0, Ti=0, steady 1D imply T_sensor=q(R+gap/k)>=qR.
    A subzero ice interface, transient storage or 2D transfer invalidates this
    particular bound. It does not identify an intrinsic contact resistance.
    """
    values=[sensor_temperature_C,sensor_absolute_bound_C,mean_heat_flux_W_m2]
    if any(type(x) not in (int,float) or not math.isfinite(x) for x in values):raise ValueError('finite numeric inputs required')
    if sensor_absolute_bound_C<0 or mean_heat_flux_W_m2<=0 or temperate_steady_assumption_explicit is not True:
        raise ValueError('positive q and explicit steady temperate assumption required')
    upper=(sensor_temperature_C+sensor_absolute_bound_C)/mean_heat_flux_W_m2
    return dict(necessary_R_upper_bound_m2K_W=upper,
                nonnegative_R_feasible=upper>=0,
                scope='necessary bound under declared steady temperate 1D assumptions only',
                intrinsic_contact_R_identified=False)
