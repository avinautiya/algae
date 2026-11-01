"""Small synthetic coupled-model development study, not a held-out benchmark."""
import copy,json
from pathlib import Path
from model import forecast

def run(root):
    c=json.loads((root/'example.json').read_text());out={}
    variants={'base':{},'warm_ice':{'ice_subcooling_K':0.},'bare_face':{'face_cover_fraction':0.},
              'nonzero_transmission':{'aperture_m':.005},
              'contact_slope_low':{'contact_friction_slope':.005},
              'contact_slope_high':{'contact_friction_slope':.02}}
    for name,updates in variants.items():
        config=copy.deepcopy(c);config['physics'].update(updates)
        try:out[name]=dict(status='complete',**forecast(config))
        except ValueError as exc:out[name]=dict(status='refused',reason=str(exc),configuration=config)
    refinement=[]
    for dt in (1800.,900.,450.):
        config=copy.deepcopy(c);config['step_s']=dt;r=forecast(config)
        refinement.append(dict(step_s=dt,face_melt_m3=r['totals']['face_melt_m3'],
                               floor_melt_m3=r['totals']['floor_melt_m3'],**r['checks']))
    result=dict(evidence='SYNTHETIC DEVELOPMENT ONLY; all variants specified after component development; no professional-skill claim',
                variants=out,timestep_refinement=refinement)
    (root/'development_results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(variants={k:dict(status=v['status'],face_melt_m3=v['totals']['face_melt_m3']) if v['status']=='complete' else v for k,v in out.items()},refinement=refinement),indent=2))
    return result

def figures(root,data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for name in ('base','warm_ice','bare_face','nonzero_transmission'):
        if data['variants'][name]['status']!='complete':continue
        rows=data['variants'][name]['series'];t=[r['time_s']/3600 for r in rows]
        axes[0,0].plot(t,[r['toe_fraction'] for r in rows],label=name)
        axes[0,1].plot(t,[r['face_melt_m3'] for r in rows],label=name)
    rows=data['variants']['base']['series'];t=[r['time_s']/3600 for r in rows]
    for k in ('floor_melt_m3','bypass_melt_m3'):axes[1,0].plot(t,[r[k] for r in rows],label=k)
    refs=data['timestep_refinement'];axes[1,1].plot([r['step_s']/60 for r in refs],[r['face_melt_m3'] for r in refs],'o-')
    for ax in axes.flat:ax.grid(alpha=.2)
    axes[0,0].set_ylabel('Toe share of prescribed inflow');axes[0,0].legend(fontsize=8)
    axes[0,1].set_ylabel('Cumulative face melt (m³)');axes[0,1].legend(fontsize=8)
    axes[1,0].set_ylabel('Cumulative channel melt (m³)');axes[1,0].legend(fontsize=8)
    for ax in (axes[0,0],axes[0,1],axes[1,0]):ax.set_xlabel('Time (hours)')
    axes[1,1].set_xlabel('Coupling step (minutes)');axes[1,1].set_ylabel('24-hour face melt (m³)')
    refused=[k for k,v in data['variants'].items() if v['status']!='complete']
    fig.suptitle('Synthetic coupled development study — not field validation\nRefused configurations: '+(', '.join(refused) or 'none'),fontsize=11);fig.tight_layout()
    fig.savefig(root/'coupled_development.png',dpi=180);fig.savefig(root/'coupled_development.pdf');plt.close(fig)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plot',action='store_true');a=p.parse_args()
    root=Path(__file__).resolve().parent;data=run(root)
    if a.plot:figures(root,data)
