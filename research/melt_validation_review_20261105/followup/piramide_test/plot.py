"""Render saved matched scores with both declared weighting schemes."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent

def main():
    r=json.loads((R/'PUBLISHED_RESULTS.json').read_text())
    names=['robust_candidate','KP1_2026_S22_fixed_effects','pdd_over_h','DETI_reconstruction_original_current_gate_TT0','DETI_reconstruction_original_current_gate_TT1','DETI_reconstruction_modified_lagged_gate_TT0','DETI_reconstruction_modified_lagged_gate_TT1']
    labels=['Robust candidate','Published KP1','Degree-day / thickness','DETI original gate, 0°C','DETI original gate, 1°C','DETI lagged gate, 0°C','DETI lagged gate, 1°C']
    data=[next(s for s in r['scores'] if not s['label_end'] and s['method']==n) for n in names]
    plt.rcParams.update({'font.size':8,'font.family':'DejaVu Sans','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(7.08,4.0),sharey=True,constrained_layout=True)
    for ax,key,title in zip(axes,['equal_site_MAE_m','pooled_day_MAE_m'],['a  Equal-site weighting','b  Equal-day weighting']):
        ax.barh(np.arange(len(data)),[s[key]*1000 for s in data],color=['#0072B2','#D55E00','#777777']+['#009E73']*4)
        ax.set_yticks(np.arange(len(data)),labels);ax.set_xlim(0,4);ax.set_xlabel('Daily MAE (mm ice)');ax.set_title(title,loc='left',weight='bold')
    axes[0].invert_yaxis()
    fig.suptitle('Pirámide 2022–2023 external campaign\n97 matched stake-days • 3 sites • observed meteorology\nDETI bars are equation reconstructions; no confidence intervals',fontsize=10)
    fig.savefig(R/'comparison.svg');fig.savefig(R/'comparison.png',dpi=220);plt.close(fig)
    with (R/'figure_source.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['method','equal_site_MAE_m','pooled_day_MAE_m','site1_days','site2_days','site3_days']);w.writeheader()
        for s in data:w.writerow(dict(method=s['method'],equal_site_MAE_m=s['equal_site_MAE_m'],pooled_day_MAE_m=s['pooled_day_MAE_m'],**{f'site{k}_days':s['sites'][str(k)]['n'] for k in [1,2,3]}))
    (R/'FIGURE_CAPTION.md').write_text('Figure. Error of frozen models on matched camera-read stake windows at Pirámide. Left: each of three sites is weighted equally. Right: each of97complete days is weighted equally. Sites contribute48,40and9days, respectively. The plotted subset also requires the preceding lagged forcing for every DETI reconstruction. All methods use exactly the same site/date windows under the start-labelled hourly-mean convention. DETI bars are explicitly equation reconstructions with0/1°Cthreshold and original/lagged threshold gates, not independently verified author-code reproductions. No statistical confidence intervals are shown. These are correlated days in one external campaign, not97independent glacier replications. The end-labelled alternative and every model remain in PUBLISHED_RESULTS.json.\n')

if __name__=='__main__':main()
