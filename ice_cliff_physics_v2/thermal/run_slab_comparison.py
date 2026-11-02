"""Bounded, unfit comparison of one-node closure and independent slab PDE."""
from pathlib import Path
import csv,json,math,hashlib
from wet_debris import ThermalState,ThermalParameters,advance as lumped_step
from resolved_slab import SlabState,SlabParameters,advance as slab_step

ROOT=Path(__file__).resolve().parent
P=SlabParameters(.05,1.,1.8e6,1.)
TAU=P.diffusion_time_s


def simulate(biot,n=64,dt_fraction=.0005,horizon_fraction=1.,initial_C=0.,flow_capacity=400.,inlet_C=.5):
    C=P.volumetric_capacity_J_m3_K*P.area_m2*P.thickness_m
    K=P.conductivity_W_m_K*P.area_m2/P.thickness_m
    G=biot*K
    state=SlabState(tuple([initial_C]*n));lumped=ThermalState(initial_C)
    lp=ThermalParameters(C,K,0.,G);t=0.;ice=lic=0.;residual=0.;samples=[]
    events=[x*TAU for x in [.01,.1,1.] if x<=horizon_fraction];stop=TAU*horizon_fraction;seen=set()
    while t<stop-1e-10:
        dt=min(dt_fraction*TAU,stop-t)
        for e in events:
            if t+1e-8<e<t+dt:dt=e-t
        r=slab_step(state,P,dt,inlet_C,flow_capacity,0.,G);state=r['state']
        l=lumped_step(lumped,lp,dt,inlet_C,flow_capacity);lumped=l.state
        ice+=r['covered_ice_heat_J'];lic+=l.covered_ice_heat_J;residual+=r['global_energy_residual_J'];t+=dt
        if any(abs(t-e)<1e-7 for e in events) and round(t/TAU,6) not in seen:
            seen.add(round(t/TAU,6))
            H=flow_capacity*(-math.expm1(-G/flow_capacity)) if flow_capacity>0 else 0.
            steady=K*H*inlet_C/(K+H) if H>0 else 0.
            scale=steady*t if steady>0 else C*initial_C
            samples.append(dict(Bi_nominal=biot,time_s=t,time_over_diffusion_time=t/TAU,
                slab_covered_heat_J=ice,lumped_covered_heat_J=lic,
                lumped_minus_slab_J=lic-ice,
                signed_difference_over_reference_heat=None if abs(ice)<1e-9 else (lic-ice)/ice,
                difference_over_steady_flux_time_or_initial_storage=None if scale==0 else (lic-ice)/scale,
                reference_grid_cells=n,reference_step_s=dt_fraction*TAU,
                slab_discrete_energy_residual_J=residual))
    return dict(samples=samples,covered_heat_J=ice,lumped_heat_J=lic,state=list(state.temperatures_C),
                accumulated_discrete_energy_residual_J=residual,final_time_s=t)


def main():
    cases=[]
    for bi in [.1,1.,10.]:cases.extend(simulate(bi)['samples'])
    hot=simulate(0.,horizon_fraction=1.,initial_C=4.,flow_capacity=0.,inlet_C=0.)
    cold_storage=P.volumetric_capacity_J_m3_K*P.area_m2*P.thickness_m*4
    for row in hot['samples']:
        t=row['time_over_diffusion_time']
        remaining=sum(8/(math.pi**2*(2*j+1)**2)*math.exp(-((2*j+1)*math.pi)**2*t/4) for j in range(300))
        row['analytic_dry_slab_heat_J']=cold_storage*(1-remaining)
        row['relative_numerical_error_to_analytic_dry_slab']=(row['slab_covered_heat_J']-row['analytic_dry_slab_heat_J'])/row['analytic_dry_slab_heat_J']
    mesh=[dict(cells=n,**{k:v for k,v in simulate(1.,n=n,dt_fraction=.00025,horizon_fraction=.2).items() if k not in ['state','samples']}) for n in [16,32,64,128]]
    time=[dict(step_fraction=dt,**{k:v for k,v in simulate(1.,n=128,dt_fraction=dt,horizon_fraction=.2).items() if k not in ['state','samples']}) for dt in [.01,.005,.0025,.00125,.000625]]
    result=dict(classification='unfit synthetic closure applicability experiment, not field accuracy or historical novelty',
        physical_parameters=P.__dict__,diffusion_time_s=TAU,
        same_lumped_capacity_and_ice_conductance='C=rhoc*A*d, K=k*A/d; not fitted to reference',
        declared_shared_water_closure='quasi-steady finite-capacity parallel exchange; slab diffusion implementation independent',
        main_cases=cases,dry_hot_initial_pulse=hot['samples'],spatial_refinement=mesh,time_refinement=time,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'wet_debris.py',ROOT/'resolved_slab.py',Path(__file__)]})
    (ROOT/'SLAB_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n')
    with (ROOT/'SLAB_COMPARISON.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(cases[0]));w.writeheader();w.writerows(cases)
    print(json.dumps(dict(main_cases=cases,dry_hot_pulse=hot['samples']),indent=2))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
        for bi in [.1,1.,10.]:
            rows=[r for r in cases if r['Bi_nominal']==bi]
            axes[0].plot([r['time_over_diffusion_time'] for r in rows],[100*r['difference_over_steady_flux_time_or_initial_storage'] for r in rows],'o-',label=f'h d/k = {bi:g}')
        axes[0].axhline(0,color='k',lw=.8);axes[0].set_xscale('log');axes[0].set(xlabel='Time / slab diffusion time',ylabel='(Lumped − resolved heat) / (steady flux × time), %',title='Same physical parameters; no fit')
        rows=hot['samples'];axes[1].plot([r['time_over_diffusion_time'] for r in rows],[r['slab_covered_heat_J']/cold_storage for r in rows],'o-',label='Resolved slab')
        axes[1].plot([r['time_over_diffusion_time'] for r in rows],[r['lumped_covered_heat_J']/cold_storage for r in rows],'s--',label='One-node closure')
        axes[1].set_xscale('log');axes[1].set(xlabel='Time / slab diffusion time',ylabel='Initial debris heat transferred to ice, fraction',title='Hot debris pulse: dry, insulated top')
        for ax in axes:ax.legend(fontsize=8);ax.grid(alpha=.25)
        fig.suptitle('Energy conservation does not establish debris transient accuracy',fontsize=12)
        for suffix in ['png','pdf']:fig.savefig(ROOT/f'slab_closure_comparison.{suffix}',dpi=180)
        plt.close(fig)
    except ImportError:
        pass

if __name__=='__main__':main()
