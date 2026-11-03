"""Reproduce all-nine-method figures from released, frozen artifacts; no fitting."""
import copy,datetime as dt,hashlib,json,os
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures'
OUT.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'.mplcache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from anomaly_index import load_site,dose

NAMES={'anomaly_exp_air':'Exponential moisture–air (primary)',
       'air_pdd':'Air temperature','buried_pdd':'Buried temperature',
       'gradient_index':'Positive-gradient index','constant_rate':'Constant rate',
       'shifted_ratio_air':'Shifted-ratio moisture–air','centered_softplus_air':'Centered-softplus moisture–air',
       'calendar_exp_air':'Calendar-adjusted air','moisture_only_exp':'Moisture-only exponential'}
COLORS={'anomaly_exp_air':'#ac3535','air_pdd':'#5b6770','calendar_exp_air':'#17658a'}
SELECT=['anomaly_exp_air','air_pdd','calendar_exp_air']
NOTES='One glacier; four serial periods per stake. Empirical conditional hindcasts, not weather forecasts.\nCamera uncertainty is unknown; no error bars or confidence intervals are inferred.'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=260,facecolor='white')
    fig.savefig(OUT/(name+'.svg'),facecolor='white')
    plt.close(fig)


def annotate_note(fig,text=NOTES):
    fig.text(.02,.012,text,fontsize=8,color='#404040',va='bottom')


def main():
    folder=ROOT/'frozen_tower2_final'
    result=json.loads((folder/'RESULTS.json').read_text());issued=json.loads((folder/'PREDICTIONS.json').read_text())
    freeze=json.loads((folder/'FREEZE.json').read_text());development=json.loads((ROOT/'DEVELOPMENT_FIT.json').read_text())
    if sha(folder/'PREDICTIONS.json')!=result['prediction_sha256']:raise ValueError('result/prediction artifact mismatch')
    metrics=result['model_metrics'];models=freeze['models']
    if set(models)!=set(metrics) or len(models)!=9:raise ValueError('must plot all nine frozen models')
    # Recompute metrics from the existing period records as a plotting check.
    for method in models:
        rows=metrics[method]['period_errors'];mae=sum(abs(x['predicted_mmwe']-x['observed_mmwe']) for x in rows)/4
        if len(rows)!=4 or abs(mae-metrics[method]['mean_absolute_period_ablation_error_mmwe'])>1e-10:
            raise ValueError('period/MAE discrepancy')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
                         'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    order=sorted(models,key=lambda x:metrics[x]['mean_absolute_period_ablation_error_mmwe'])
    fig,ax=plt.subplots(figsize=(9.7,6.4));fig.subplots_adjust(left=.34,right=.94,bottom=.18,top=.85)
    values=[metrics[m]['mean_absolute_period_ablation_error_mmwe'] for m in order]
    ax.barh(range(9),values,color=[COLORS.get(m,'#b6c0c6') for m in order],height=.70)
    ax.set_yticks(range(9),[NAMES[m] for m in order]);ax.invert_yaxis()
    for i,value in enumerate(values):ax.text(value+.35,i,f'{value:.2f}',va='center',fontsize=9)
    ax.set_xlim(0,max(values)*1.14);ax.set_xlabel('Mean absolute period ablation error (mm water equivalent)')
    ax.set_title('Tower2 transfer test — all nine frozen methods',loc='left',pad=25)
    ax.text(0,1.025,'Lower is better; the prespecified primary is highlighted in red.',transform=ax.transAxes,fontsize=9,color='#404040')
    ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True);annotate_note(fig)
    save(fig,'figure1_all_nine_mae')

    period_labels=['16–30 Dec\n2022','30 Dec–13 Jan\n2022–23','13–27 Jan\n2023','27 Jan–10 Feb\n2023']
    reference=metrics['anomaly_exp_air']['period_errors'];observed=np.array([r['observed_mmwe'] for r in reference]);x=np.arange(4)
    fig,axs=plt.subplots(2,1,figsize=(9.2,7),sharex=True,gridspec_kw={'height_ratios':[1.6,1.]})
    fig.subplots_adjust(left=.12,right=.97,top=.90,bottom=.21,hspace=.14)
    axs[0].plot(x,observed,'ko-',lw=2.2,label='Observed stake ablation',zorder=5)
    for method in SELECT:
        values=np.array([r['predicted_mmwe'] for r in metrics[method]['period_errors']])
        axs[0].plot(x,values,'s--',color=COLORS[method],lw=1.6,label=NAMES[method])
        axs[1].plot(x,values-observed,'o-',color=COLORS[method],lw=1.5)
    axs[0].set_ylabel('Period ablation\n(mm water equivalent)');axs[0].legend(fontsize=8,ncol=2,loc='upper left')
    axs[0].set_ylim(0,max(max(observed),max(r['predicted_mmwe'] for r in reference))*1.30)
    axs[0].set_title('Tower2 predictions and observed ablation',loc='left')
    axs[1].axhline(0,color='black',lw=.9);axs[1].set_ylabel('Prediction − observation\n(mm water equivalent)')
    axs[1].set_xticks(x,period_labels);axs[1].set_xlabel('Exact source-local evaluation periods; shared 03:00 endpoints')
    for ax in axs:ax.grid(axis='y',alpha=.18)
    annotate_note(fig);save(fig,'figure2_period_predictions')

    # Fixed parameters, fixed known local TRAIN amplitude: recomputation only,
    # never an optimization or a new parameter fit.
    feature1=load_site('Tower1');y0=development['train_mmwe'];ydev=np.array(development['development_mmwe'])
    dev={};dev_predictions={}
    for method in models:
        exposure=dose(feature1,method,freeze['parameters'][method]);pred=y0*exposure[1:]/exposure[0]
        dev[method]=float(np.mean(np.abs(pred-ydev)));dev_predictions[method]=pred.tolist()
        if method in development['fits'] and abs(dev[method]-development['fits'][method]['selected']['development_MAE_mmwe'])>1e-9:
            raise ValueError('fixed development recomputation mismatch')
    fig,ax=plt.subplots(figsize=(9.7,6.8));fig.subplots_adjust(left=.34,right=.94,bottom=.22,top=.79)
    y=np.arange(9);height=.33
    ax.barh(y-height/2,[dev[m] for m in order],height,color='#a4b4bf',label='Tower1 development (used to learn global parameters)')
    ax.barh(y+height/2,[metrics[m]['mean_absolute_period_ablation_error_mmwe'] for m in order],height,
            color=[COLORS.get(m,'#5c788c') for m in order],label='Tower2 later periods (parameters frozen)')
    ax.set_yticks(y,[NAMES[m] for m in order]);ax.invert_yaxis();ax.set_xlabel('Mean absolute period ablation error (mm water equivalent)')
    ax.set_title('Development and transfer errors for all nine methods',loc='left',pad=57)
    ax.legend(loc='lower left',bbox_to_anchor=(0,1.005),fontsize=8,frameon=False);ax.grid(axis='x',alpha=.18);ax.set_axisbelow(True)
    annotate_note(fig,'Separate local TRAIN amplitude per stake; global parameters learned only on Tower1.\nOne glacier; four serial periods per stake; conditional hindcasts, not weather forecasts.\nCamera uncertainty is unknown; no error bars or confidence intervals are inferred.')
    save(fig,'figure3_development_transfer')

    feature2=load_site('Tower2');beta=freeze['parameters']['anomaly_exp_air']
    fig,axs=plt.subplots(1,2,figsize=(11,4.9));fig.subplots_adjust(left=.10,right=.98,bottom=.25,top=.85,wspace=.30)
    site_colors={'Tower1':'#7d5c99','Tower2':'#ba6030'};period_factors={}
    for f in (feature1,feature2):
        label=f['tower'];days={}
        for stamp,value in zip(f['labels'][:-1],f['theta']-f['theta_reference']):days.setdefault(stamp.date(),[]).append(float(value))
        dates=sorted(days);daily=[np.mean(days[t]) for t in dates]
        axs[0].plot(dates,daily,color=site_colors[label],lw=1.5,label=label+(' (development)' if label=='Tower1' else ' (transfer)'))
        expdose=dose(f,'anomaly_exp_air',beta);plain=dose(f,'air_pdd',0.)
        factors=(expdose[1:]/expdose[0])/(plain[1:]/plain[0]);period_factors[label]=factors.tolist()
        axs[1].plot(x,factors,'o-',color=site_colors[label],label=label)
    axs[0].axhline(0,c='black',lw=.8);axs[0].set_ylabel('SM3 minus own TRAIN mean\n(reported volumetric-fraction units)')
    axs[0].set_title('A  Observed moisture anomalies',loc='left');axs[0].legend(fontsize=8)
    import matplotlib.dates as mdates
    axs[0].xaxis.set_major_locator(mdates.WeekdayLocator(interval=2));axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    axs[0].set_xlabel('December 2022–February 2023 (source-local dates)')
    axs[1].axhline(1,c='black',lw=.8);axs[1].set_ylabel('Exponential-index prediction / air-index prediction')
    axs[1].set_xticks(x,['Period1','Period2','Period3','Period4']);axs[1].set_title('B  Effect of the frozen moisture term',loc='left');axs[1].legend(fontsize=8)
    for ax in axs:ax.grid(axis='y',alpha=.15)
    annotate_note(fig,'Daily means are descriptive, not uncertainty intervals. Same fixed β=5.4; no refitting.\nRaw SM3 is a sensor proxy with flagged source limits; trends do not prove a causal moisture effect.\nTower1 uses 30 Dec 00:00; Tower2 uses 30 Dec 03:00 at their shared period boundary.')
    save(fig,'figure4_moisture_transfer')

    numbers=dict(model_order=order,development_MAE_mmwe=dev,
        Tower2_MAE_mmwe={m:metrics[m]['mean_absolute_period_ablation_error_mmwe'] for m in models},
        primary='anomaly_exp_air',primary_unchanged=True,observed_Tower2_period_mmwe=observed.tolist(),
        development_predictions_fixed_parameters=dev_predictions,moisture_to_plain_prediction_factors=period_factors,
        field_parameters_refitted=False,all_nine_models_plotted=True,
        camera_uncertainty='unknown; no error bars/confidence intervals',
        scope='one glacier, two stakes, four serial evaluation periods; empirical conditional hindcast')
    (OUT/'FIGURE_NUMBERS.json').write_text(json.dumps(numbers,indent=2)+'\n')
    inputs=[folder/'FREEZE.json',folder/'PREDICTIONS.json',folder/'RESULTS.json',ROOT/'DEVELOPMENT_FIT.json']
    manifest=dict(script_sha256=sha(Path(__file__)),source_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
        figures=[p.name for p in sorted(OUT.glob('figure*.png'))],fitting_performed=False,visual_inspection_pending=True)
    (OUT/'FIGURE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':main()
