"""Reproducible plots of locked field comparison; no model fitting."""
import json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/glacier_round2_plot_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    folder=ROOT/'model/frozen_piramide_v3/melt_comparison'
    data=json.loads((folder/'RESULTS.json').read_text());p=json.loads((folder/'PREDICTIONS.json').read_text())
    output=ROOT/'figures';output.mkdir(exist_ok=True)
    names=[data['primary_candidate']]+p['baseline_ids']
    labels=['Forward thermal\n(candidate)','Derivative thermal','Constant rate','Buried-temperature\nindex','Air-temperature index\n(gaps reconstructed)']
    obs=data['model_metrics'][names[0]]['observed_total_mmwe']
    vals=[data['model_metrics'][n]['predicted_total_mmwe'] for n in names]
    plt.rcParams.update({'font.size':9,'savefig.dpi':160})
    fig,ax=plt.subplots(figsize=(9,4.8))
    ax.scatter(range(len(names)),vals,c=['#146d78']+['#6f7b83']*4,s=45)
    ax.axhline(obs,color='black',label=f'Native camera target: {obs:.2f} mm w.e.')
    ax.set(xticks=range(5),xticklabels=labels,ylabel='Total ablation (mm water equivalent)',title='Pirámide Tower 1 — locked later-period comparison')
    ax.legend();fig.subplots_adjust(left=.11,right=.97,bottom=.26,top=.90)
    fig.text(.02,.02,'One tower, one coverage-amended evaluation block; camera precision is unknown.\nNo confidence intervals. Empirical training gains are separate from the raw physical energy ledger.',fontsize=8)
    fig.savefig(output/'field_predictions.png');fig.savefig(output/'field_predictions.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4.8))
    errs=[data['model_metrics'][n]['absolute_total_ablation_error_mmwe'] for n in names]
    ax.set_ylim(0,max(errs)*1.17)
    bars=ax.bar(range(5),errs,color=['#146d78']+['#6f7b83']*4)
    ax.bar_label(bars,fmt='%.2f',padding=3);ax.set(xticks=range(5),xticklabels=labels,ylabel='Absolute cumulative error (mm w.e.)',title='Predeclared primary metric — all baseline comparisons retained')
    fig.subplots_adjust(left=.11,right=.97,bottom=.26,top=.90)
    fig.text(.02,.02,'Nominal errors against the same native endpoint operator; not uncertainty-qualified superiority.\nThe derivative method is OLS, not a reproduction of the source NTLS/full energy-balance model.',fontsize=8)
    fig.savefig(output/'field_absolute_errors.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.5,4.8))
    colors={'dry':'#b17336','proxy':'#146d78','saturated':'#577cba'}
    for method,marker,x in [('forward','o',0),('derivative_OLS','s',1)]:
        for capacity,color in colors.items():
            rows=[v['absolute_total_ablation_error_mmwe'] for n,v in data['model_metrics'].items() if n.startswith(method+'__') and n.endswith('__'+capacity)]
            xs=[x+(i-(len(rows)-1)/2)*.025 for i in range(len(rows))]
            ax.scatter(xs,rows,color=color,marker=marker,label=capacity if x==0 else None)
    ax.set(xticks=[0,1],xticklabels=['Forward thermal','Derivative thermal'],ylabel='Absolute cumulative error (mm w.e.)',title='Declared parameter scenarios — not independent observations')
    ax.legend(title='Capacity assumption');fig.subplots_adjust(bottom=.22)
    fig.text(.02,.02,'Each point is one geometry × forcing reconstruction × capacity scenario on the SAME field event.\nThe primary remains middle-gap / left reconstruction / proxy capacity; no best-case selection.',fontsize=8)
    fig.savefig(output/'field_scenarios.png');plt.close(fig)
if __name__=='__main__':main()
