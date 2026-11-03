"""Reproducible scientific figure from Tower1 TRAIN-development diagnostics."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent


def main():
    data=json.loads((ROOT/'TRAINING_DIAGNOSIS.json').read_text())
    extra=json.loads((ROOT/'ADDITIONAL_IDENTIFIABILITY.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(11.5,8),layout='constrained')
    grid=sorted({r['k_upper_W_mK'] for r in data['all_candidates']})
    heat=np.array([[next(r['unweighted_temperature_MSE_C2'] for r in data['all_candidates'] if r['k_upper_W_mK']==u and r['k_lower_W_mK']==l) for u in grid] for l in grid])
    a=axes[0,0];im=a.imshow(np.log10(heat),origin='lower',aspect='auto',cmap='viridis')
    a.set_xticks(range(len(grid)),grid);a.set_yticks(range(len(grid)),grid)
    a.set_xlabel('Upper-layer conductivity (W m⁻¹ K⁻¹)');a.set_ylabel('Lower-layer conductivity (W m⁻¹ K⁻¹)')
    a.set_title('A  Two-layer TRAIN temperature loss',loc='left')
    win=next(w for w in data['winners'] if w['family']=='two_layer' and w['objective']=='unweighted_temperature_MSE_C2')
    a.plot(grid.index(win['k_upper_W_mK']),grid.index(win['k_lower_W_mK']),'x',color='white',ms=10,mew=2)
    fig.colorbar(im,ax=a,label='log₁₀ MSE (°C²)',shrink=.85)
    a=axes[0,1]
    selected=[('Homogeneous MSE',next(w for w in data['winners'] if w['family']=='homogeneous' and w['objective']=='unweighted_temperature_MSE_C2')),
              ('Two-layer MSE',win),('Near-basal objective',next(w for w in data['winners'] if w['family']=='two_layer' and w['objective']=='basal_temperature_and_gradient_loss'))]
    for name,w in selected:a.plot(range(4),w['sensor_RMSE_C'],'o-',label=name)
    a.axhline(.4,c='gray',ls='--',lw=1,label='0.4°C sensor reference')
    a.set_xticks(range(4),['DT2','DT3','DT4','DT5']);a.set_ylabel('TRAIN RMSE (°C)')
    a.set_title('B  Better temperature fit still has discrepancy',loc='left');a.legend(fontsize=8)
    a=axes[1,0]
    scale=[r for r in extra['known_ice_absolute_scale_tests'] if r['family']=='two_layer']
    a.bar(['Dry','Probe proxy','Saturated'],[r['raw_TRAIN_melt_mmwe'] for r in scale],color=['#7b99b2','#346b98','#163e63'])
    target=win['already_inspected_TRAIN_camera_mmwe'];a.axhline(target,c='#b64335',ls='--',label=f'Known TRAIN ablation: {target:.1f} mmwe')
    a.set_ylabel('Unscaled TRAIN melt (mm water equivalent)')
    a.set_title('C  Same temperatures, different absolute heat',loc='left')
    a.text(.04,.95,'Temperature change < 3 × 10⁻¹³ °C',va='top',transform=a.transAxes,fontsize=9)
    a.legend(fontsize=8,loc='upper right',bbox_to_anchor=(1,.80))
    a=axes[1,1]
    runs=extra['finite_difference_step_sensitivity'];h=np.array([r['log_parameter_difference_step'] for r in runs]);s=np.array([r['singular_values'][-1] for r in runs])
    a.loglog(h,s,'o-',label='Smallest computed singular value')
    a.loglog(h,s[0]*(h/h[0])**2,'--',c='gray',label='h² differentiation error')
    a.set_xlabel('Log-parameter difference step h');a.set_ylabel('Scaled singular value')
    a.set_title('D  Apparent scale information is numerical',loc='left');a.legend(fontsize=8)
    fig.suptitle('Tower1 TRAIN development — no new field validation',fontsize=14)
    output=ROOT/'figures';output.mkdir(exist_ok=True)
    fig.savefig(output/'training_physics_diagnosis.png',dpi=180)
    fig.savefig(output/'training_physics_diagnosis.svg')
    plt.close(fig)


if __name__=='__main__':main()
