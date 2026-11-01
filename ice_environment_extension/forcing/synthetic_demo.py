"""Synthetic warm/cold bare/debris comparison; no observational inputs.

Run: python3 synthetic_demo.py --output synthetic_ledger.json [--figure plot.png]
Matplotlib is needed only when --figure is requested.
"""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
from energy import (BareIceState, ColumnState, ColumnParameters,
                    step_bare_ice, step_column)


def simulate():
    latent = 334000.0
    parameters = ColumnParameters(100000.0, 5.0, latent)
    states = {'warm_bare': BareIceState(0), 'cold_bare': BareIceState(2e6),
              'warm_debris': ColumnState(0, 0), 'cold_debris': ColumnState(0, 2e6)}
    initial = {name: asdict(state) for name, state in states.items()}
    records = []
    totals = {name: {k: 0.0 for k in ('atmospheric_j_m2', 'water_j_m2',
              'debris_storage_j_m2', 'cold_storage_j_m2', 'latent_j_m2',
              'residual_j_m2')} for name in states}
    # Three synthetic days: square-wave daytime warming and nighttime cooling.
    for hour in range(72):
        atmosphere = 180.0 if 6 <= hour % 24 < 18 else -60.0
        for name, state in states.items():
            if isinstance(state, BareIceState):
                result = step_bare_ice(state, latent, 3600, atmosphere, 10.0)
            else:
                result = step_column(state, parameters, 3600, atmosphere, 10.0)
            states[name] = result.state
            for key in totals[name]:
                totals[name][key] += getattr(result, key)
            records.append({'hour_end': hour + 1, 'scenario': name,
                            'atmospheric_w_m2': atmosphere, 'water_w_m2': 10.0,
                            **asdict(result)})
    return {'label': 'SYNTHETIC illustrative fixed-area scenarios; not site predictions',
            'assumptions': 'Departed meltwater; lumped cold reserve; zero-C debris interface; unlimited ice.',
            'parameters': asdict(parameters), 'initial_states': initial,
            'final_states': {k: asdict(v) for k, v in states.items()},
            'integrated_ledgers_j_m2': totals, 'records': records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='synthetic_ledger.json')
    parser.add_argument('--figure')
    args = parser.parse_args()
    data = simulate()
    Path(args.output).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    if args.figure:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
        for name in data['final_states']:
            rows = [r for r in data['records'] if r['scenario'] == name]
            axes[0].plot([r['hour_end'] for r in rows],
                         [r['state']['melted_kg_m2'] for r in rows], label=name)
            axes[1].plot([r['hour_end'] for r in rows],
                         [r['state']['cold_reserve_j_m2'] / 1e6 for r in rows])
        axes[0].set_ylabel('Cumulative melt (kg/m²)')
        axes[1].set_ylabel('Cold reserve (MJ/m²)')
        axes[1].set_xlabel('Synthetic time (hours)')
        axes[0].legend()
        fig.suptitle('SYNTHETIC fixed-area forcing comparison')
        fig.tight_layout()
        fig.savefig(args.figure, dpi=150)
        plt.close(fig)
    print(args.output)


if __name__ == '__main__':
    main()
