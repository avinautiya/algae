"""Synthetic conditional model demonstration, not a glacier result."""
import copy
import json
from pathlib import Path
from .model import forecast


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    route=dict(ice_z_m=0.,width_m=.5,ice_area_m2=1.,heat_transfer_W_m2_K=200.,
               debris_thickness_m=0.,debris_conductivity_W_m_K=1.)
    c=dict(pond_area_m2=10.,pond_floor_m=-1.,initial_stage_m=.189,
           initial_temperature_K=.5,horizon_s=15*86400.,output_step_s=6*3600.,
           routes=[copy.deepcopy(route),copy.deepcopy(route)],
           forcing=[dict(time_s=0.,inflow_m3_s=.07,temperature_K=.5)],
           evidence='SYNTHETIC reduced model; prescribed clearing, not predicted grain transport')
    c['routes'][0]['debris_thickness_m']=.2
    out=Path(__file__).parent/'records';out.mkdir(exist_ok=True)
    (Path(__file__).parent/'sample_input.json').write_text(json.dumps(c,indent=2)+'\n')
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    summary={}
    for name,day in [('no clearing',None),('clear day 1',1),('clear day 10',10)]:
        case=copy.deepcopy(c)
        if day is not None:case['removals']=[dict(time_s=day*86400.,route=0)]
        result=forecast(case)
        (out/(name.replace(' ','_')+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        time=[r['time_s']/86400 for r in result['series']]
        axs[0].plot(time,[r['toe_fraction'] for r in result['series']],label=name)
        axs[1].plot(time,[r['ice_z_m'][0]-r['ice_z_m'][1] for r in result['series']],label=name)
        axs[2].plot(time,[r['stage_m'] for r in result['series']],label=name)
        summary[name]=dict(final_toe_fraction=result['series'][-1]['toe_fraction'],
                           ice_latent_energy_J=result['ledger']['ice_latent_energy_J'],
                           rhs_evaluations=result['rhs_evaluations'],checks=result['checks'])
    for ax in axs:ax.set_xlabel('Time (days)')
    axs[0].set_ylabel('Toe-route water share')
    axs[1].set_ylabel('Toe minus bypass ice elevation (m)')
    axs[2].set_ylabel('Reservoir water level (m)')
    axs[0].legend(fontsize=8)
    fig.suptitle('Synthetic two-control forecast — prescribed clearing; no spatial or field validation',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(out/f'dynamic_demonstration.{ext}',dpi=180)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
