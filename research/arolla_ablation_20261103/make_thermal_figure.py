"""Plot declared synthetic example only; never reads field targets."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
def main():
    data=json.loads((ROOT/'model/thermal_memory/example_results.json').read_text())
    fig,axes=plt.subplots(2,1,figsize=(8,6.5),sharex=True)
    for i,r in enumerate(data['scenarios']):
        history=r['history'];t=[h['time_s']/3600 for h in history]
        axes[0].plot(t,[h['q_basal_W_m2'] for h in history],label=f'Declared initial profile {i+1}')
        total=0.;m=[]
        for h in history:
            total+=h['step_s']*h['q_basal_W_m2']/r['latent_heat_J_kg'];m.append(total)
        axes[1].plot(t,m)
    axes[0].set(ylabel='Basal heat delivery (W/m²)',title='Synthetic example of stored debris heat affecting later melt')
    axes[0].legend();axes[1].set(xlabel='Simulated time (hours)',ylabel='Cumulative conductive melt (mm w.e.)')
    fig.subplots_adjust(left=.13,right=.97,top=.91,bottom=.17,hspace=.1)
    fig.text(.02,.015,'Identical imposed surface forcing and thermal parameters; only initial debris temperature differs.\nNumerical demonstration, not field validation. Two scenarios do not bound all environmental uncertainty.',fontsize=8)
    (ROOT/'figures').mkdir(exist_ok=True);fig.savefig(ROOT/'figures/thermal_memory_synthetic.png',dpi=170)
if __name__=='__main__':main()
