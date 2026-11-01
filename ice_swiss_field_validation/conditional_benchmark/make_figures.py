"""Paper figures from saved conditional results; no fitting or selection."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main(source,out):
    result=json.loads(source.read_text());nom=result['scenarios']['clock+0h_T+0.00C_D1.30m'];out.mkdir(parents=True,exist_ok=True)
    methods=['constant_rate_fitted','constant_h_fitted','db04_unfitted','db03_unfitted','gnielinski_smooth_unfitted','db04_scale_fitted']
    labels=['Fitted constant rate','Fitted constant h','Unfitted DB0.4','Unfitted DB0.3','Unfitted smooth Gnielinski','Fitted scaled DB0.4']
    fig,axes=plt.subplots(1,2,figsize=(11,4.7),constrained_layout=True)
    for ax,site in zip(axes,['P3','P5']):
        values=[nom['test_metrics'][site][m]['MAE_m'] for m in methods]
        ax.barh(labels,values,color=['#aaa','#8e8e8e','#348abd','#58b7d0','#467821','#d95f02'])
        ax.invert_yaxis();ax.set_xlabel('Later-interval MAE (m)');ax.set_title(f'{site}: n={nom["test_metrics"][site][methods[0]]["n"]} dependent intervals');ax.set_xlim(0,1.65);ax.grid(axis='x',alpha=.2)
    fig.suptitle('Swiss channel component comparison under imposed daily-mean forcing\nRetrospective temporal test — one stream, no general superiority established',fontsize=12)
    for ext in ['png','pdf']:fig.savefig(out/f'conditional_nominal_scores.{ext}',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),constrained_layout=True)
    for ax,site in zip(axes,['P3','P5']):
        for base,label,color in [('constant_h_fitted','Versus fitted constant h','#348abd'),('gnielinski_smooth_unfitted','Versus unfitted smooth Gnielinski','#467821')]:
            xs=[];ys=[]
            for clock in [-2,0,2]:
                for shift in [-.05,0,.05]:
                    s=result['scenarios'][f'clock{clock:+d}h_T{shift:+.2f}C_D1.30m'];m=s['test_metrics'][site]
                    xs.append(len(xs));ys.append(m[base]['MAE_m']-m['db04_scale_fitted']['MAE_m'])
            ax.plot(xs,ys,'o-',label=label,color=color)
        ax.axhline(0,color='black',lw=1);ax.set_title(site);ax.set_ylabel('Baseline MAE − fitted scaled DB MAE (m)');ax.set_xlabel('Fixed clock/temperature scenario index');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Engineering sensitivity at D=1.30 m\nPositive means scaled DB improves MAE; scenarios are not independent replicates',fontsize=11)
    for ext in ['png','pdf']:fig.savefig(out/f'conditional_sensitivity.{ext}',dpi=180)
    plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.source,a.out)
