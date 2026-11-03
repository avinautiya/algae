"""Render only frozen retrospective results. Never refits a model."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures'; OUT.mkdir(exist_ok=True)
def main():
    d=json.loads((ROOT/'model/results/frozen_ablation_comparison/RESULTS.json').read_text())
    names=['constant_rate','positive_degree_day','net_shortwave_dose','published_DEBCF_component_unresolved']
    labels=['Constant rate','Degree day\n(predeclared primary)','Net shortwave','DEBCF*\n(unresolved)']
    s=d['model_summary']
    plt.rcParams.update({'font.size':10,'savefig.dpi':170})
    fig,ax=plt.subplots(figsize=(9,4.7),layout='constrained')
    for i,n in enumerate(names):
        a=s[n]; p=a['primary_nominal_prediction_mmwe']; lo,hi=a['primary_timing_scenario_prediction_range_mmwe']
        ax.errorbar(i,p,yerr=[[p-lo],[hi-p]],fmt='o',capsize=6,color='#176878' if i<3 else '#777777')
    ax.axhline(315,color='black',label='Measured ice ablation: 315 mm w.e.')
    ax.set(xticks=range(4),xticklabels=labels,ylabel='Ice ablation (mm water equivalent)',title='Arolla H1D retrospective temporal test')
    ax.legend(loc='upper left')
    fig.text(.02,.01,'Bars span 81 visit-time/bin scenarios; not confidence intervals. Source measurement uncertainty is additional.\nDEBCF output is contextual: ice/snow component unresolved. One glacier and one evaluation weather interval.',fontsize=8)
    ax.set_position([.10,.26,.85,.61])
    fig.savefig(OUT/'arolla_primary.png');fig.savefig(OUT/'arolla_primary.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.6),layout='constrained')
    for j,n in enumerate(names[:3]):
        rows=[r for r in d['rows'] if r['model']==n and r['nominal']]
        ax.plot([r['stake'] for r in rows],[r['signed_error_mmwe'] for r in rows],'-o',label=labels[j].split('\n')[0])
    ax.axhline(0,color='black',lw=.8);ax.set(ylabel='Prediction minus observation (mm w.e.)',title='All four stakes share the same weather episode')
    ax.legend();fig.savefig(OUT/'arolla_signed_errors.png');plt.close(fig)
    audit=json.loads((ROOT/'timing/FIELD_RANKING_AUDIT.json').read_text())
    fig,ax=plt.subplots(figsize=(8,4.6),layout='constrained')
    interpretations=list(audit['status_counts'])
    counts=[sum(audit['status_counts'][x].values()) for x in interpretations]
    bars=ax.bar(['If ±0.1 m per reading','If ±0.1 m per interval'],counts,color='#c9823a')
    ax.bar_label(bars);ax.set(ylim=(0,280),ylabel='Unresolved pair/scenario rankings',title='No robust pairwise ranking under either bounded interpretation')
    fig.text(.02,.01,'243 = 3 model pairs × 81 timing scenarios, not independent observations.\nSource uncertainty association is unspecified; neither interpretation is a statistical confidence interval.',fontsize=8)
    ax.set_position([.13,.23,.80,.63]);fig.savefig(OUT/'arolla_ranking_uncertainty.png');plt.close(fig)
if __name__=='__main__':main()
