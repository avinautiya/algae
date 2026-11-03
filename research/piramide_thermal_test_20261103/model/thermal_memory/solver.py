"""One-dimensional conductive debris over melting ice; no refreezing.

Coordinate z points downward. q_surface enters the debris and q_basal leaves it.
Cell-centered FV with backward Euler, exact Dirichlet half-cell resistances.
Source declarations are gates, not automatic certification of their truth.
"""
from dataclasses import dataclass
import math

INDEPENDENT = {'independently_measured', 'independent_literature'}
SYNTHETIC = 'declared_synthetic'


def positive(value, name, *, zero=False):
    if type(value) not in (float, int) or not math.isfinite(value):
        raise ValueError(name + ' must be finite numeric (not bool)')
    if value < 0 or (value == 0 and not zero):
        raise ValueError(name + (' must be nonnegative' if zero else ' must be positive'))
    return float(value)


def provenance(record, name, synthetic):
    if not isinstance(record, dict):
        raise ValueError(name + ': provenance object required')
    allowed = INDEPENDENT | ({SYNTHETIC} if synthetic else set())
    if record.get('status') not in allowed:
        raise ValueError(name + ': target-fitted, unknown or unsupported provenance refused')
    for key in ('source', 'location', 'method'):
        if not isinstance(record.get(key), str) or not record[key].strip():
            raise ValueError(name + ': nonempty ' + key + ' required')
    if record.get('uses_evaluation_melt') is not False:
        raise ValueError(name + ': uses_evaluation_melt must explicitly be false')


@dataclass(frozen=True)
class Column:
    conductivity_W_mK: float
    volumetric_heat_capacity_J_m3K: float
    depth_m: float
    initial_temperature_C: tuple
    sources: dict


@dataclass(frozen=True)
class SurfaceBoundary:
    """Nonnegative *debris-surface* temperatures, native-interval values held piecewise.

    value[i] applies to [i*interval_seconds, (i+1)*interval_seconds).
    Air temperature alone is not this boundary condition.
    """
    surface_temperature_C: tuple
    source: dict
    interval_seconds: float = 3600.


def _tridiagonal(lower, diagonal, upper, rhs):
    """Thomas elimination for the positive-definite conduction operator."""
    n = len(diagonal)
    b, r = list(diagonal), list(rhs)
    for j in range(1, n):
        multiplier = lower[j-1] / b[j-1]
        b[j] -= multiplier * upper[j-1]
        r[j] -= multiplier * r[j-1]
    result = [0.] * n
    result[-1] = r[-1] / b[-1]
    for j in range(n-2, -1, -1):
        result[j] = (r[j] - upper[j] * result[j+1]) / b[j]
    if not all(math.isfinite(x) for x in result):
        raise ValueError('nonfinite thermal solve')
    return result


def simulate(column, boundary, *, horizon_s, max_horizon_s, dt_s,
             synthetic=False, latent_heat_J_kg=333700., record_profiles=False):
    """Return per-unit-plan-area heat and its latent-equivalent ice ablation.

    The ice remains at 0 C, with no basal conductive loss, refreezing or moving
    geometry. Melt mm w.e. = basal J/m2 / latent J/kg, since 1 kg/m2=1 mm w.e.
    This predicts a sub-debris conductive component, not full cliff melt.
    """
    if type(record_profiles) is not bool:
        raise ValueError('record_profiles must be bool')
    interval = positive(boundary.interval_seconds, 'native boundary interval')
    if type(synthetic) is not bool:
        raise ValueError('synthetic must be bool')
    k = positive(column.conductivity_W_mK, 'conductivity')
    C = positive(column.volumetric_heat_capacity_J_m3K, 'heat capacity')
    depth = positive(column.depth_m, 'depth')
    horizon = positive(horizon_s, 'horizon')
    cap = positive(max_horizon_s, 'maximum permitted horizon')
    step_limit = positive(dt_s, 'dt')
    latent = positive(latent_heat_J_kg, 'latent heat')
    if horizon > cap:
        raise ValueError('horizon exceeds explicitly permitted maximum')
    if not isinstance(column.initial_temperature_C, (tuple, list)) or not column.initial_temperature_C:
        raise ValueError('explicit nonempty cell-center initial profile required')
    T = [positive(t, 'initial temperature (refreezing unsupported)', zero=True)
         for t in column.initial_temperature_C]
    if not isinstance(column.sources, dict):
        raise ValueError('column source records required')
    for name in ('conductivity', 'heat_capacity', 'depth', 'initial_profile', 'basal_ice_state'):
        provenance(column.sources.get(name), name, synthetic)
    if column.sources['heat_capacity'].get('quantity') != 'bulk_debris_volumetric_heat_capacity':
        raise ValueError('bulk debris heat capacity required, not rock heat capacity alone')
    basal_state = column.sources['basal_ice_state']
    if basal_state.get('state') != 'temperate_melting_interface' or type(basal_state.get('temperature_C')) not in (int, float) or basal_state['temperature_C'] != 0.:
        raise ValueError('source-supported temperate melting ice at 0 C required')
    provenance(boundary.source, 'surface temperature', synthetic)
    if boundary.source.get('quantity') != 'debris_surface_temperature':
        raise ValueError('boundary must be debris_surface_temperature, not air temperature')
    if not isinstance(boundary.surface_temperature_C, (tuple, list)) or not boundary.surface_temperature_C:
        raise ValueError('explicit contiguous native surface boundary required')
    surface = [positive(t, 'surface temperature (refreezing unsupported)', zero=True)
               for t in boundary.surface_temperature_C]
    if horizon > len(surface)*interval:
        raise ValueError('surface boundary does not cover requested horizon')
    # Reject numerically unrepresentable/nonphysical coefficients instead of
    # allowing an uninformative singular solve or a zero-progress time loop.
    n = len(T)
    dz = depth/n
    storage = C*dz
    interior = k/dz
    boundary_g = 2*interior
    if not all(math.isfinite(x) and x > 0 for x in (dz, storage, interior, boundary_g)):
        raise ValueError('unrepresentable discretization coefficients')
    initial_storage = storage*math.fsum(T)
    heat_in, heat_out, cumulative_residual = [], [], []
    t = 0.
    history = []
    while t < horizon:
        hour = min(int(t//interval), len(surface)-1)
        h = min(step_limit, horizon-t, (hour+1)*interval-t)
        if h <= 0 or t+h == t:
            raise ValueError('dt too small to advance time')
        capacity_rate = storage/h
        if not math.isfinite(capacity_rate):
            raise ValueError('dt creates unrepresentable capacity coefficient')
        lower = [-interior]*(n-1)
        upper = [-interior]*(n-1)
        diagonal = [capacity_rate + (interior if j else boundary_g)
                    + (interior if j < n-1 else boundary_g) for j in range(n)]
        rhs = [capacity_rate*v for v in T]
        rhs[0] += boundary_g*surface[hour]
        if not all(math.isfinite(x) for x in diagonal+rhs):
            raise ValueError('unrepresentable matrix or boundary forcing')
        new = _tridiagonal(lower, diagonal, upper, rhs)
        if min(new) < 0:
            raise ValueError('negative solution in no-refreezing scope')
        q_in = boundary_g*(surface[hour]-new[0])
        q_out = boundary_g*new[-1]
        energy_in = h*q_in
        energy_out = h*q_out
        change = storage*math.fsum(b-a for a, b in zip(T, new))
        residual = change - energy_in + energy_out
        if not all(math.isfinite(x) for x in (energy_in, energy_out, change, residual)):
            raise ValueError('unrepresentable energy ledger')
        heat_in.append(energy_in)
        heat_out.append(energy_out)
        cumulative_residual.append(residual)
        t += h
        T = new
        history.append(dict(time_s=t, step_s=h, surface_temperature_C=surface[hour],
                            q_surface_W_m2=q_in, q_basal_W_m2=q_out,
                            storage_change_J_m2=change,
                            energy_residual_J_m2=residual))
        if record_profiles:
            history[-1]['temperature_C'] = list(T)
    basal = math.fsum(heat_out)
    supplied = math.fsum(heat_in)
    final_storage = storage*math.fsum(T)
    if not all(math.isfinite(x) for x in (initial_storage, final_storage, basal, supplied)):
        raise ValueError('unrepresentable total energy ledger')
    global_residual = final_storage - initial_storage - supplied + basal
    scale = max(abs(initial_storage), abs(final_storage), abs(supplied), abs(basal), 1.)
    return dict(scope='fixed-thickness 1D sub-debris conductive ice-ablation component',
                synthetic=synthetic,
                column_parameters=dict(conductivity_W_mK=k,
                                       volumetric_heat_capacity_J_m3K=C,
                                       depth_m=depth, cells=n,
                                       basal_temperature_C=0.),
                time_options=dict(horizon_s=horizon,max_horizon_s=cap,dt_s=step_limit),
                surface_boundary_temperature_C=surface,
                surface_interval_seconds=interval,
                initial_temperature_C=list(column.initial_temperature_C),
                final_temperature_C=T, cell_center_depth_m=[(i+.5)*dz for i in range(n)],
                basal_heat_J_m2=basal, melt_mm_water_equivalent=basal/latent,
                latent_heat_J_kg=latent, surface_heat_J_m2=supplied,
                initial_sensible_storage_J_m2=initial_storage,
                final_sensible_storage_J_m2=final_storage,
                discrete_energy_residual_J_m2=global_residual,
                discrete_energy_relative_residual=global_residual/scale,
                sum_step_residual_J_m2=math.fsum(cumulative_residual),
                max_absolute_step_residual_J_m2=max(abs(x) for x in cumulative_residual),
                history=history, declared_sources=column.sources,
                surface_source=boundary.source,
                validated_against_field_melt=False, mathematical_novelty_claim=False)


def ensemble(scenarios, boundary, **options):
    """Finite declared-scenario envelope, not a continuous uncertainty proof."""
    if not isinstance(scenarios, (tuple, list)) or not scenarios:
        raise ValueError('nonempty explicit scenarios required')
    results = [simulate(column, boundary, **options) for column in scenarios]
    values = [x['melt_mm_water_equivalent'] for x in results]
    return dict(scenarios=results, declared_scenario_count=len(results),
                melt_mm_water_equivalent_range=[min(values), max(values)],
                envelope_scope='only these finite parameter/initial-state scenarios',
                continuous_parameter_bound_established=False,
                statistical_confidence_interval=False)
