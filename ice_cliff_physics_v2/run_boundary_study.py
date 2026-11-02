"""Frozen fixed-network comparison against an independent Robin-slab solution."""
import argparse,hashlib,json,math,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from thermal.resolved_slab import SlabParameters,SlabState,advance
from thermal.wet_debris import ThermalParameters,ThermalState,advance as rc_advance
from thermal.robin_reference import RobinSlab
from mechanism.wet_resistance import prediction
from run_study import atomic,digest,source_files


def finite_volume(case,cells,steps_per_tau):
    p=SlabParameters(case['thickness_m'],case['conductivity_W_m_K'],case['volumetric_heat_capacity_J_m3_K'],case['wet_area_m2'])
    tau=p.diffusion_time_s;step=tau/steps_per_tau;state=SlabState((0.,)*cells)
    totals=dict(ice=0.,water=0.,air=0.);samples={};maximum_residual=0.
    for index in range(1,10*steps_per_tau+1):
        r=advance(state,p,step,case['inlet_superheat_K'],case['water_capacity_W_K'],
                  debris_water_conductance_W_K=case['film_conductance_W_K'],atmospheric_power_W=case['absorbed_source_W'])
        state=r['state'];totals['ice']+=r['covered_ice_heat_J'];totals['water']+=r['water_to_debris_heat_J'];totals['air']+=r['atmospheric_heat_J']
        maximum_residual=max(maximum_residual,abs(r['global_energy_residual_J']),abs(r['water_exchange_residual_J']))
        if index in [steps_per_tau//10,steps_per_tau,10*steps_per_tau]:
            samples[str(round(index/steps_per_tau,3))]=dict(time_s=index*step,ice_heat_J=totals['ice'],
                ice_power_W=r['covered_ice_heat_J']/step,water_heat_J=totals['water'],atmospheric_heat_J=totals['air'],
                profile_C=state.temperatures_C)
    return dict(cells=cells,steps_per_diffusion_time=steps_per_tau,samples=samples,max_discrete_energy_residual_J=maximum_residual)

def run(out,commit):
    freeze=json.loads((ROOT/'study/FREEZE.json').read_text())
    if freeze['sources']!=source_files() or freeze['wet_boundary_cases_sha256']!=digest(ROOT/'mechanism/wet_boundary_cases.json'):
        raise ValueError('source or boundary cases changed after freeze')
    if len(commit)!=40 or any(c not in '0123456789abcdef' for c in commit):raise ValueError('published freeze commit required')
    cases=json.loads((ROOT/'mechanism/wet_boundary_cases.json').read_text());index=[]
    for c in cases:
        start=time.perf_counter();path=out/(c['id']+'.json')
        if path.exists():raise ValueError('retain old boundary outputs; use a new run directory')
        ca=c['volumetric_heat_capacity_J_m3_K']*c['thickness_m'];K=c['ice_conductance_W_K'];C=ca*c['wet_area_m2'];tau=C/K
        ref=RobinSlab.from_contact(C,K,c['water_capacity_W_K'],c['film_conductance_W_K'],c['inlet_superheat_K'],c['absorbed_source_W'],terms=512)
        ref256=RobinSlab.from_contact(C,K,c['water_capacity_W_K'],c['film_conductance_W_K'],c['inlet_superheat_K'],c['absorbed_source_W'],terms=256)
        analytic={};rc={}
        for factor in [.1,1.,10.]:
            key=str(factor);analytic[key]=ref.response(tau*factor)
            analytic[key]['series_256_vs_512_ice_difference_J']=ref256.response(tau*factor)['ice_heat_J']-analytic[key]['ice_heat_J']
            r=rc_advance(ThermalState(0.),ThermalParameters(C,K,0.,c['film_conductance_W_K']),tau*factor,
                c['inlet_superheat_K'],c['water_capacity_W_K'],atmospheric_power_W=c['absorbed_source_W'])
            rc[key]=dict(ice_heat_J=r.covered_ice_heat_J,storage_J=r.debris_storage_change_J,
                error_J=r.covered_ice_heat_J-analytic[key]['ice_heat_J'])
        fv=[finite_volume(c,n,steps) for n,steps in [(16,200),(32,200),(32,400),(64,400)]]
        steady=prediction(water_capacity_W_K=c['water_capacity_W_K'],film_conductance_W_K=c['film_conductance_W_K'],
            ice_conductance_W_K=K,inlet_superheat_K=c['inlet_superheat_K'],absorbed_external_heat_W=c['absorbed_source_W'])
        record=dict(id=c['id'],case=c,freeze_commit=commit,scope='controlled thermal-network verification, not glacier observations',
                    analytic_robin=analytic,rc=rc,finite_volume=fv,steady_reference=steady,runtime_s=time.perf_counter()-start)
        atomic(path,record);index.append(dict(id=c['id'],runtime_s=record['runtime_s']))
        print(c['id'],'completed',round(record['runtime_s'],3),flush=True)
    atomic(out/'RUN_INDEX.json',dict(freeze_commit=commit,completed=len(index),cases=index))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--freeze-commit',required=True);a=p.parse_args();run(a.out,a.freeze_commit)
