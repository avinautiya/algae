"""Synthetic prescribed-interface-flux verification, never a field prediction."""
from dataclasses import asdict
import argparse
import json
import math
from pathlib import Path
from face_profile import Profile, advance, RHO_ICE, LATENT

ROOT=Path(__file__).resolve().parent
WET_FLUX=2000.
AIR_FLUX=60.
FLOOR_FLUX=300.
WATER_STAGE=.5
HORIZON=86400.


def run_case(n_bands,dt):
    g=Profile.uniform(0.,10.,-2.,2.,0.,2.,1.,n_bands)
    initial=g;t=0.;rows=[];maximum_residual=0.;total_face=total_floor=0.
    while t<HORIZON:
        step=min(dt,HORIZON-t);energies=[]
        for b in g.bands:
            wet=max(0.,min(b.z_high_m,WATER_STAGE)-b.z_low_m)
            air=b.z_high_m-b.z_low_m-wet
            energies.append((WET_FLUX*wet+AIR_FLUX*air)*step)
        floor_energy=FLOOR_FLUX*g.floor_area_m2*step
        g,ledger=advance(g,energies,floor_energy)
        maximum_residual=max(maximum_residual,abs(ledger['volume_residual_m3']))
        total_face+=ledger['face_melt_m3'];total_floor+=ledger['floor_melt_m3'];t+=step
        rows.append(dict(time_s=t,geometry=asdict(g),floor_area_m2=g.floor_area_m2,
                         total_melt_m3=total_face+total_floor,**g.undercut_diagnostic()))
    return dict(n_bands=n_bands,dt_s=dt,initial=asdict(initial),final=asdict(g),series=rows,
                total_melt_m3=total_face+total_floor,face_melt_m3=total_face,floor_melt_m3=total_floor,
                remaining_volume_loss_m3=initial.ice_volume_m3-g.ice_volume_m3,
                max_step_volume_residual_m3=maximum_residual,
                final_polygon=g.polygon(),final_floor_z_m=g.floor_z_m,
                final_floor_area_m2=g.floor_area_m2,**g.undercut_diagnostic())


def analytic():
    aw=WET_FLUX/(RHO_ICE*LATENT);aa=AIR_FLUX/(RHO_ICE*LATENT);b=FLOOR_FLUX/(RHO_ICE*LATENT);t=HORIZON
    return dict(wet_front_x_m=2+aw*t,air_front_x_m=2+aa*t,floor_z_m=-b*t,
                total_melt_m3=(.5*aw+1.5*aa+2*b)*t+aw*b*t*t,
                interpretation='exact geometry for prescribed piecewise constant face speed and constant floor speed; not a thermal PDE solution')


def analytic_uniform():
    aw=WET_FLUX/(RHO_ICE*LATENT);aa=AIR_FLUX/(RHO_ICE*LATENT);b=FLOOR_FLUX/(RHO_ICE*LATENT);t=HORIZON
    dx=aw*t+(aa-aw)*1.5/b*math.log1p(b*t/2.)
    return dict(front_x_m=2+dx,floor_z_m=-b*t,total_melt_m3=2*dx+b*t*(2+dx),
                interpretation='exact solution of the reduced uniform-recession model, not the nonuniform physical profile')


def build():
    cases={'uniform_face':run_case(1,900.),'height_resolved_face':run_case(64,900.)}
    refinements=[run_case(64,dt) for dt in [3600.,1800.,900.]]
    exact=analytic()
    for case in refinements:
        case['floor_error_m']=abs(case['final_floor_z_m']-exact['floor_z_m'])
        case['melt_error_m3']=abs(case['total_melt_m3']-exact['total_melt_m3'])
        del case['series']
    return dict(evidence='synthetic prescribed-latent-flux geometry verification; no field accuracy or fracture prediction',
                forcing=dict(wet_face_latent_flux_W_m2=WET_FLUX,air_face_latent_flux_W_m2=AIR_FLUX,
                             floor_latent_flux_W_m2=FLOOR_FLUX,horizontal_ledge_latent_flux_W_m2=0.,
                             water_stage_m=WATER_STAGE,horizon_s=HORIZON,
                             note='externally prescribed NET LATENT fluxes, not a solved finite water heat or atmospheric energy budget'),
                cases=cases,analytic_reference=exact,analytic_uniform_model=analytic_uniform(),time_refinement=refinements)


def plots(data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(13,4.8),constrained_layout=True)
    initial=[(2,0),(2,2)]
    axes[0].plot([x for x,z in initial],[z for x,z in initial],':',color='black',label='Initial face')
    for key,label,color in [('uniform_face','Uniform recession','#d97706'),('height_resolved_face','Height-resolved recession','#007a87')]:
        c=data['cases'][key];g=c['final'];bands=g['bands']
        points=[]
        for b in reversed(bands):
            points += [(b['front_x_m'],b['z_high_m']),(b['front_x_m'],b['z_low_m'])]
        points += [(0.,bands[0]['z_low_m'])]
        axes[0].plot([x for x,z in points],[z for x,z in points],color=color,label=label)
        axes[1].plot([r['time_s']/3600 for r in c['series']],[r['floor_area_m2'] for r in c['series']],color=color,label=label)
    axes[0].axhline(WATER_STAGE,color='#436dba',lw=.7,ls='--',label='Prescribed stage')
    axes[0].set(xlim=(1.8,2.8),ylim=(-.15,2.05),xlabel='Front position x (m)',ylabel='Elevation z (m)',title='Nonuniform forcing permits an undercut')
    axes[0].legend(fontsize=7)
    axes[1].set(xlabel='Elapsed time (h)',ylabel='Basal footprint (m²)',title='Same initial body and imposed flux law')
    ref=data['time_refinement'];axes[2].loglog([r['dt_s']/3600 for r in ref],[r['melt_error_m3'] for r in ref],'o-')
    axes[2].set(xlabel='Coupling step (h)',ylabel='Volume error vs analytic geometry (m³)',title='Resolved profile approaches exact solution')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Prescribed-flux geometry benchmark — no flow, collapse or field validation',fontsize=12)
    for suffix in ['png','pdf']:fig.savefig(ROOT/f'profile_benchmark.{suffix}',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--plot',action='store_true');a=p.parse_args()
    d=build();(ROOT/'EXAMPLE_RESULTS.json').write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
    if a.plot:plots(d)
    print(json.dumps(dict(analytic=d['analytic_reference'],cases={k:{n:v for n,v in c.items() if n in ['total_melt_m3','final_floor_area_m2','maximum_horizontal_undercut_m']} for k,c in d['cases'].items()},refinement=[{k:r[k] for k in ['dt_s','floor_error_m','melt_error_m3']} for r in d['time_refinement']]),indent=2))
