"""Plot saved field scores without rerunning, fitting, or selecting a model."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent

def main():
    a = json.loads((ROOT/'event_benchmark/FIELD_RESULTS.json').read_text())
    b = json.loads((ROOT/'event_benchmark/FIELD_RESULTS_2023.json').read_text())
    years = [(2022, a['metrics']), (2023, b['scenarios']['base']['metrics'])]
    keys = list(a['metrics']['pooled']['models'])
    labels = ['2021 constant', '2021 degree-day', 'Fixed energy balance', '2021 fitted energy balance']
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    for ax, (year, metrics) in zip(axes, years):
        sites = list(metrics)
        x = np.arange(len(sites))
        for j, (key, label) in enumerate(zip(keys, labels)):
            vals = [metrics[s]['models'][key]['MAE_h'] for s in sites]
            ax.bar(x+(j-1.5)*.19, vals, width=.19, label=label)
        ax.set_xticks(x, [s.replace('P_IC_', '') for s in sites])
        ax.set_ylabel('Event-duration MAE (hours)')
        ax.set_title(f'{year} observed-weather hindcast')
        ax.set_ylim(0, 14)
        ax.grid(axis='y', alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle('Public Kennicott stake data: skill varies by year and cliff\nTime to next 10 cm marker; atmospheric components only', fontsize=12)
    out = ROOT/'figures'; out.mkdir(exist_ok=True)
    for suffix in ['png', 'pdf']:
        fig.savefig(out/f'field_scores.{suffix}', dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    for key, label in zip(keys, labels):
        values = [b['scenarios'][s]['metrics']['pooled']['models'][key]['MAE_h'] for s in ['SW_minus15percent', 'base', 'SW_plus15percent']]
        ax.plot([-15, 0, 15], values, 'o-', label=label)
    ax.set(xlabel='Engineering change in shortwave forcing (%)', ylabel='Pooled event-duration MAE (hours)', title='2023 forcing sensitivity — scenarios, not confidence bounds')
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    for suffix in ['png', 'pdf']:
        fig.savefig(out/f'forcing_sensitivity.{suffix}', dpi=180)
    plt.close(fig)

if __name__ == '__main__':
    main()
