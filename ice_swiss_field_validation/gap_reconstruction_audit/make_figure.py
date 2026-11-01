"""Plot all frozen masking results, without method selection or model fitting."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    d=json.loads(a.results.read_text());a.out.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    labels={'linear_endpoint_bridge':'Linear endpoint bridge','constant_previous_endpoint':'Previous endpoint','constant_next_endpoint':'Next endpoint'}
    for ax,station in zip(axes,['P3','P5']):
        for method,label in labels.items():
            ax.plot([1,6,24],[d['scores'][f'{station}_{h}h'][method]['MAE_C'] for h in [1,6,24]],'o-',label=label)
        ax.set(xscale='log',xticks=[1,6,24],xticklabels=['1','6','24'],xlabel='Withheld block duration (hours)',ylabel='MAE of block-mean temperature (°C)',title=station)
        ax.grid(alpha=.2)
        ax.text(.02,.97,'Blocks: '+', '.join(str(d['counts'][f'{station}_{h}h']['retained_blocks']) for h in [1,6,24]),transform=ax.transAxes,va='top',fontsize=8)
    axes[0].legend(fontsize=8,loc='lower right')
    fig.suptitle('Controlled temperature-gap reconstruction\nNative-data masking; correlated blocks; not melt-model validation',fontsize=11)
    for suffix in ['png','pdf']:fig.savefig(a.out/f'gap_reconstruction.{suffix}',dpi=180)
    plt.close(fig)

if __name__=='__main__':main()
