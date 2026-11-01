"""Conditional finite-flood exclusion certificate, not a sufficient rescue rule.

Uses ideal broad-crested controls, a cube sliding criterion and the uniform
thermal resistance of the model. It does not apply to lift, rolling, granular
avalanching, spatially patchy insulation or mechanical ice collapse.
"""
import math

def cube_critical_depth(D, friction, drag, density_ratio):
    if not all(math.isfinite(x) for x in (D,friction,drag,density_ratio)) or D<=0 or friction<0 or drag<=0 or density_ratio<=1:
        raise ValueError('invalid cube inputs')
    x=(-friction+math.sqrt(friction**2+2*drag*friction*density_ratio))/drag
    if x>1:x=2*friction*(density_ratio-1)/drag
    return D*x


def minimum_sliding_discharge(c,gap,D):
    dc=cube_critical_depth(D,c.friction,c.drag,c.grain_density/c.water_density)
    head=1.5*dc
    K=(2/3)**1.5*math.sqrt(c.gravity)*c.width
    return K*(head**1.5+max(gap+head,0.)**1.5)


def exclude_rescue(c, *, toe_ice_z, bypass_ice_z, protected_solid_stock_m3,
                   grain_D_m, peak_discharge_m3_s, duration_s):
    vals=(toe_ice_z,bypass_ice_z,protected_solid_stock_m3,grain_D_m,peak_discharge_m3_s,duration_s)
    if not all(math.isfinite(x) for x in vals) or min(protected_solid_stock_m3,grain_D_m,peak_discharge_m3_s,duration_s)<=0:
        raise ValueError('invalid flood state')
    As=(1-c.porosity)*c.width*c.length
    coarse_depth=protected_solid_stock_m3/As
    # Bare-water conductance is a conservative ceiling even for patchy cover.
    h_eff=c.h_water
    max_lowering=h_eff*c.water_temperature*duration_s/(c.ice_density*c.latent)
    gap_lower=toe_ice_z+coarse_depth-bypass_ice_z-max_lowering
    required=minimum_sliding_discharge(c,gap_lower,grain_D_m)
    return dict(coarse_mobilization_excluded=peak_discharge_m3_s<required,
                optimistic_minimum_sliding_discharge_m3_s=required,
                most_favorable_gap_lower_bound_m=gap_lower,
                optimistic_toe_lowering_m=max_lowering,
                warning='conditional exclusion only; failure to exclude is not proof of rescue')


def thermal_lockout(c, *, toe_ice_z, bypass_ice_z, peak_discharge_m3_s,
                    duration_s, return_discharge_m3_s=None,
                    extra_toe_heat_W_m2=0.):
    """Even perfect instantaneous clearing cannot restore returned ordinary flow.

    Optimistically allocate heat to the toe and omit all bypass lowering. This
    only certifies the declared fixed-width vertical-incision route geometry.
    extra_toe_heat is a nonnegative upper bound on other delivered ice heat.
    """
    Qreturn=c.discharge if return_discharge_m3_s is None else return_discharge_m3_s
    vals=(toe_ice_z,bypass_ice_z,peak_discharge_m3_s,duration_s,Qreturn,extra_toe_heat_W_m2)
    if not all(math.isfinite(x) for x in vals) or min(peak_discharge_m3_s,Qreturn)<=0 or min(duration_s,extra_toe_heat_W_m2)<0:
        raise ValueError('invalid lockout inputs')
    K=(2/3)**1.5*math.sqrt(c.gravity)*c.width
    # A generous heat bound: assign all upstream sensible heat to toe; also
    # enforce the bare-water Robin flux ceiling. Additional heating is explicit.
    advective=c.water_density*c.water_cp*peak_discharge_m3_s*c.water_temperature/(c.width*c.length)
    water_flux=min(advective,c.h_water*c.water_temperature)
    max_lowering=(water_flux+extra_toe_heat_W_m2)*duration_s/(c.ice_density*c.latent)
    min_gap=toe_ice_z-bypass_ice_z-max_lowering
    dry_head=(Qreturn/K)**(2/3)
    return dict(reconnection_excluded_even_after_perfect_clearing=min_gap>=dry_head,
                optimistic_remaining_ice_gap_m=min_gap,ordinary_flow_activation_head_m=dry_head,
                optimistic_maximum_toe_lowering_m=max_lowering,
                scope='conditional finite-event exclusion; geometric bypass creation, collapse and widening absent')
