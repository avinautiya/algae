"""Locked retrospective field results only: no fitting or new target access."""
import json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/glacier_round2_plot_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
    future=json.loads((ROOT/'model/frozen_piramide_future_batch_final/melt_comparison/RESULTS.json').read_text())
    prior=json.loads((ROOT/'model/frozen_piramide_v3/melt_comparison/RESULTS.json').read_text())
    names=[future['primary_candidate']]+list(future['primary_comparisons'])
    labels=['Forward thermal\n(candidate)','Derivative thermal','Constant rate','Buried-temperature\nindex','Air-temperature index\n(gaps reconstructed)']
    out=ROOT/'figures';out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':9,'savefig.dpi':170})
    values=[future['model_metrics'][n]['mean_absolute_period_ablation_error_mmwe']for n in names]
    fig,ax=plt.subplots(figsize=(9,4.7));bars=ax.bar(range(5),values,color=['#146d78']+['#7d8790']*4)
    ax.bar_label(bars,fmt='%.2f',padding=3);ax.set(ylim=(0,max(values)*1.2),xticks=range(5),xticklabels=labels,ylabel='Mean absolute period error (mm w.e.)',title='Three later periods — unchanged parameters and training gains')
    fig.subplots_adjust(left=.11,right=.97,bottom=.27,top=.90)
    fig.text(.02,.015,'Three serial periods at one tower, not independent glaciers. No confidence intervals.\nAll 39 models retained; fourth proposed period excluded before outcomes for missing air-temperature support.',fontsize=8)
    fig.savefig(out/'future_batch_mae.png');fig.savefig(out/'future_batch_mae.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.5,4.8))
    obs=[prior['model_metrics'][names[0]]['observed_total_mmwe']]+[r['observed_mmwe']for r in future['model_metrics'][names[0]]['period_errors']]
    ax.plot(range(4),obs,'-ok',label='Native camera observations')
    for name,color,label in [(names[0],'#146d78','Forward thermal'),(names[-1],'#b07b34','Air-temperature index')]:
        vals=[prior['model_metrics'][name]['predicted_total_mmwe']]+[r['predicted_mmwe']for r in future['model_metrics'][name]['period_errors']]
        ax.plot(range(4),vals,'-o',color=color,label=label)
    ax.axvline(.5,color='#999999',ls='--',lw=1)
    ax.set(xticks=range(4),xticklabels=['Dec16–Dec30','Dec30–Jan13','Jan13–Jan27','Jan27–Feb10'],ylabel='Period ablation (mm w.e.)',title='Earlier result retained beside the fixed later-period batch')
    ax.legend();fig.subplots_adjust(bottom=.24)
    fig.text(.02,.02,'2022–23 source-local labels; exact endpoint timestamps are in the result files.\nObserved-temperature-conditioned hindcasts; source camera precision remains unknown.',fontsize=8)
    fig.savefig(out/'four_period_predictions.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.5,4.8));colors={'dry':'#b17336','proxy':'#146d78','saturated':'#577cba'}
    for method,x,marker in [('forward',0,'o'),('derivative_OLS',1,'s')]:
        for capacity,color in colors.items():
            vv=[v['mean_absolute_period_ablation_error_mmwe']for n,v in future['model_metrics'].items()if n.startswith(method+'__')and n.endswith('__'+capacity)]
            ax.scatter([x+(i-(len(vv)-1)/2)*.025 for i in range(len(vv))],vv,c=color,marker=marker,label=capacity if x==0 else None)
    ax.axhline(values[-1],color='black',ls='--',label='Air-index primary MAE')
    ax.set(xticks=[0,1],xticklabels=['Forward thermal','Derivative thermal'],ylabel='Three-period MAE (mm w.e.)',title='Finite parameter scenarios — no post-scoring model selection')
    ax.legend();fig.subplots_adjust(bottom=.23)
    fig.text(.02,.02,'Each point uses the same three field periods; points are not independent samples.\nCapacity/geometry/reconstruction scenarios do not bound every environmental uncertainty.',fontsize=8)
    fig.savefig(out/'future_batch_scenarios.png');plt.close(fig)
if __name__=='__main__':main()
