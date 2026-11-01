"""Reproduce a labelled synthetic geometric comparison and figure."""
import json
from pathlib import Path
import numpy as np
from .bounds import opening_bound, exclusion


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(__file__).parent/'records';out.mkdir(exist_ok=True)
    # Two resolved 2 m x 1 m ice-bed perturbations, identical added ice volume.
    # These are bed elevations AFTER perfect clearing, not debris resistance.
    strip=np.zeros((40,20));strip[15:25,:]=.08
    corridor=np.zeros_like(strip);corridor[15:25,10:]=.16
    levels=np.linspace(0,.10,101)
    curves={}
    for name,z in [('transverse ice ridge',strip),('open ice corridor',corridor)]:
        curves[name]=[opening_bound(z,level_upper_m=float(h),dx_m=.05,dy_m=.05,
                     minimum_cross_section_width_m=.20)['latent_energy_lower_bound_J'] for h in levels]
    selected=opening_bound(strip,level_upper_m=.02,dx_m=.05,dy_m=.05,
                           minimum_cross_section_width_m=.20)
    record=dict(evidence='SYNTHETIC GEOMETRY, NOT A FLOW OR GLACIER RESULT',
                added_ice_volume_m3={'transverse':float(strip.sum()*.05**2),
                                    'corridor':float(corridor.sum()*.05**2)},
                selected_bound=selected,exclusion=exclusion(selected,720000.),
                level_upper_m=levels.tolist(),energy_lower_bound_J=curves)
    (out/'synthetic_bounds.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    for ax,z,title in zip(axs[:2],(strip,corridor),('Transverse ice ridge','Open ice corridor')):
        im=ax.imshow(z,origin='lower',extent=(0,1,0,2),vmin=0,vmax=.16,aspect='auto')
        ax.set(xlabel='Transverse distance (m)',ylabel='Streamwise distance (m)',title=title)
    fig.colorbar(im,ax=list(axs[:2]),label='Post-clearance ice-bed height (m)')
    for name,vals in curves.items():axs[2].plot(levels,np.array(vals)/1e6,label=name)
    axs[2].axhline(.72,color='black',ls=':',label='Illustrative heat budget')
    axs[2].set(xlabel='Upper possible water level (m)',ylabel='Required latent energy lower bound (MJ)',title='Necessary opening condition')
    axs[2].legend(fontsize=8)
    fig.suptitle('Synthetic geometry only — exclusion is not a reconnection forecast',fontsize=11)
    for suffix in ('png','pdf'):fig.savefig(out/f'synthetic_bounds.{suffix}',dpi=180)
    print(json.dumps(record['exclusion'],indent=2))


if __name__=='__main__':main()
