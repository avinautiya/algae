"""Paired field evidence accounting; no parameter fitting or model selection.

Bounds are supplied engineering/measurement intervals, never inferred confidence
intervals. A robust result is conditional on those intervals covering reality.
"""
import math


def interval(value, name):
    if len(value) != 2:
        raise ValueError(name + ' needs two endpoints')
    lo, hi = value
    if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) for v in value) or lo > hi:
        raise ValueError('invalid ' + name)
    return float(lo), float(hi)


def worst_case_advantage(observation, candidate, baseline):
    """Exact lower bound on baseline absolute error minus candidate error.

    Minimize over the observation and both independent prediction intervals.
    For fixed truth y the minimizing baseline error is distance to its interval;
    the maximizing candidate error occurs at a candidate endpoint. The remaining
    continuous piecewise-linear function changes slope only at baseline endpoints
    and the candidate midpoint, so enumeration of these and observation endpoints
    is exact. Shared systematic errors can make this conservative; no independence
    distribution is assumed. Positive means advantage for every admitted value.
    """
    yl, yh = interval(observation, 'observation')
    cl, ch = interval(candidate, 'candidate')
    bl, bh = interval(baseline, 'baseline')
    knots = [yl, yh] + [x for x in [bl, bh, (cl + ch)/2] if yl <= x <= yh]
    return min(max(bl-y, 0., y-bh) - max(abs(cl-y), abs(ch-y)) for y in knots)


def compare(records, candidate, baselines):
    """Compare all predeclared baselines within each observed history/target.

    Each record has id, history, target_type, unit, observation bounds and
    prediction bounds. No pooled mixing of distinct physical targets is allowed.
    Records sharing a history are not counted as independent replications.
    """
    if not records or not baselines or candidate in baselines or len(set(baselines)) != len(baselines):
        raise ValueError('nonempty records and distinct baseline names required')
    ids = set(); groups = {}
    for row in records:
        if row['id'] in ids:
            raise ValueError('duplicate measurement id')
        ids.add(row['id'])
        key = (row['history'], row['target_type'], row['unit'])
        if not all(isinstance(v, str) and v for v in key):
            raise ValueError('explicit history, physical target and unit required')
        interval(row['observation'], 'observation')
        for model in [candidate] + list(baselines):
            interval(row['predictions'][model], model)
        groups.setdefault(key, []).append(row)
    results = []
    for (history, target, unit), rows in sorted(groups.items()):
        comparisons = {}
        for baseline in baselines:
            advantages = [worst_case_advantage(r['observation'], r['predictions'][candidate], r['predictions'][baseline]) for r in rows]
            comparisons[baseline] = dict(n_observations=len(rows),
                lower_bound_mean_absolute_error_advantage=math.fsum(advantages)/len(rows),
                robustly_better=math.fsum(advantages)>0)
        results.append(dict(history=history, target_type=target, unit=unit, baselines=comparisons))
    return dict(candidate=candidate, groups=results,
        consistent_within_supplied_bounds=all(v['robustly_better'] for g in results for v in g['baselines'].values()),
        interpretation='conditional bound comparison; not a confidence interval, independent-site test, or full-system validation')
