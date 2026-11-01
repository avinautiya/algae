"""Reproducible synthetic checks of the reachability bound, not glacier data."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from recovery_budget import static_bound


def main():
    base = Path(__file__).parent
    p = dict(area_m2=10., coefficient=.7, exponent=1.5,
             bypass_m=0., toe_m=.2, initial_stage_m=.05)
    qc = p['coefficient']*p['toe_m']**p['exponent']
    rows = []
    for eps in (.001, .003, .01, .03, .1, .3, 1., 3., 10.):
        b = static_bound(**p, peak_m3_s=qc*(1+eps))
        rows.append(dict(excess_peak_fraction=eps, minimum_time_s=b.minimum_time_s,
                         minimum_volume_m3=b.minimum_volume_m3))
    # Piecewise rectangular pulse at the same peak; subcritical duration must fail.
    b = static_bound(**p, peak_m3_s=.1)
    checks = []
    for ratio in (.5, .9, 1.1, 2.):
        duration = ratio*b.minimum_time_s
        for erosion_m_s in (0., .0001, .001):
            def rhs(t, y):
                zb = p['bypass_m']-erosion_m_s*t
                return [(.1-p['coefficient']*max(y[0]-zb, 0.)**1.5)/10.]
            def event(t, y): return y[0]-.2
            event.terminal = True
            event.direction = 1
            sol = solve_ivp(rhs, (0, duration), [.05], events=event,
                            rtol=1e-10, atol=1e-12)
            checks.append(dict(duration_ratio=ratio, bypass_erosion_m_s=erosion_m_s,
                               reached=bool(len(sol.t_events[0])),
                               final_stage_m=float(sol.y[0,-1])))
    record = dict(evidence='synthetic declared scalar model only', parameters=p,
                  near_threshold=rows, independent_ode_checks=checks)
    (base/'synthetic_checks.json').write_text(json.dumps(record, indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.7))
    eps = np.array([r['excess_peak_fraction'] for r in rows])
    ax[0].semilogx(eps, [r['minimum_time_s'] for r in rows], 'o-')
    ax[0].set(xlabel='Peak excess above static activation threshold',
              ylabel='Minimum pulse duration (s)', title='Storage penalty near activation')
    ax[1].semilogx(eps, [r['minimum_volume_m3'] for r in rows], 'o-')
    ax[1].axhline(10*(.2-.05), color='k', ls='--', label='Storage-only volume')
    ax[1].set(xlabel='Peak excess above static activation threshold',
              ylabel='Minimum total inflow (m³)', title='Bypass leakage costs additional water')
    ax[1].legend(fontsize=8)
    fig.suptitle('Synthetic first-contact bounds — no environmental calibration', fontsize=10)
    fig.tight_layout()
    for ext in ('png', 'pdf'): fig.savefig(base/f'recovery_budget.{ext}', dpi=180)


if __name__ == '__main__': main()
