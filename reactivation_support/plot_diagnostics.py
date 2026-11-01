from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .connectivity import minimax_crest

if __name__=='__main__':
    root=Path(__file__).parent
    rows=json.loads((root/'records/smooth_control_diagnostic.json').read_text())+json.loads((root/'records/smooth_control_refinement.json').read_text())
    fig,ax=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    choices=[('abrupt, dx=.05',False,.05),('graded, dx=.05',True,.05),('graded, dx=.025',True,.025)]
    for off,(label,smooth,dx) in enumerate(choices):
        errors=[]
        for layout in ('bare','strip'):
            r=next(r for r in rows if r['layout']==layout and r['smooth_inlet']==smooth and r['dx']==dx and not r['smooth_deposit'])
            errors.append(r['numerical_face_fraction']-r['ideal_reference_fraction'])
        ax[0].bar(np.arange(2)+(off-1)*.24,errors,width=.22,label=label)
    ax[0].axhline(.05,color='black',ls='--',lw=.8);ax[0].axhline(-.05,color='black',ls='--',lw=.8)
    ax[0].axhline(0,color='gray',lw=.8)
    ax[0].set_xticks([0,1],['Bare','Strip'])
    ax[0].set(ylabel='Numerical fraction − ideal control fraction',title='Entrance isolation: different test geometries')
    ax[0].legend(fontsize=8)
    z=np.array([[0,1,1],[0,0,1],[1,0,0],[1,1,0]],float);r=minimax_crest(z)
    image=ax[1].imshow(z,origin='upper',cmap='Blues',vmin=0,vmax=1)
    path=np.array(r['path']);ax[1].plot(path[:,1],path[:,0],'o-',color='crimson',lw=2)
    ax[1].set(title='Synthetic winding-path diagnostic',xlabel='Cross-stream cell',ylabel='Streamwise cell')
    ax[1].text(.5,-.19,'Straight-column crest 1; connected-path crest 0',ha='center',transform=ax[1].transAxes,fontsize=9)
    fig.colorbar(image,ax=ax[1],label='Synthetic elevation (arbitrary units)')
    fig.suptitle('Independent development diagnostics — not final evaluation',fontsize=12)
    out=root/'figures';out.mkdir(exist_ok=True)
    fig.savefig(out/'independent_diagnostics.png',dpi=180);fig.savefig(out/'independent_diagnostics.pdf');plt.close(fig)
