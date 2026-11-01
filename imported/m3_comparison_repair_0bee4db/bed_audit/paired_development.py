"""Bounded one-worker coarse diagnostic, NEVER a frozen campaign replacement."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
from .operator import lower_ice

ROOT=Path(__file__).resolve().parents[1]


def worker(mode,out):
    sys.path.insert(0,str(ROOT))
    from snapshot import cases, swe2d
    mutations=[]
    original=swe2d.Sim.lower_bed
    def instrumented(sim,dt):
        before=sim.zb.copy()
        if mode=='legacy':
            result=original(sim,dt)
        else:
            plan=lower_ice(before,sim.ice,sim.last_qacc,dt,sim.dx**2,
                           rho_i=swe2d.RHO_I,latent_J_kg=swe2d.LF)
            sim.zb=plan['bed']
            dz=plan['lowering_m'];rm=np.zeros_like(sim.ice);rm[sim.i0:sim.i1+1]=True
            sim.last_reach=tuple(float((dz*(lane&rm)).sum()*sim.dx**2) for lane in (sim.toe,sim.byp))
            result=(plan['ice_volume_m3'],float((dz*sim.toe).sum()*sim.dx**2),float((dz*sim.byp).sum()*sim.dx**2))
        actual=float(np.where(sim.ice,before-sim.zb,0.).sum()*sim.dx**2)
        mutations.append(dict(time_day=sim.t/86400,declared_melt_m3=float(result[0]),actual_ice_lowering_m3=actual,
                              nonice_cells_changed=int(np.sum((sim.zb!=before)&(~sim.ice)))))
        return result
    swe2d.Sim.lower_bed=instrumented
    start=time.monotonic()
    result=cases.run('strip',dx=.1,days=3.5,dtm=21600.,hyd_max=60.)
    result['snapshot_maxrss_raw_labelled_MB']=result['maxrss_MB']
    if sys.platform=='darwin':result['maxrss_MB']/=1024.
    result['memory_units_note']='raw pinned ru_maxrss/1024 retained; Darwin value normalized to MiB'
    lost=cases.contact_times(result['records'])[0]
    record=dict(evidence='EXPLORATORY coarse paired operator sensitivity, not production data',mode=mode,
                model_change='ice-only unbounded lowering; removing the artificial ice floor requires boundary/domain revalidation' if mode!='legacy' else 'pinned legacy lower_bed',
                wall_s=time.monotonic()-start,loss_day=lost,operator_updates=mutations,
                declared_minus_actual_melt_m3=sum(x['declared_melt_m3']-x['actual_ice_lowering_m3'] for x in mutations),
                run=result)
    out.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--worker',choices=['legacy','ice_only'])
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.worker:
        worker(a.worker,a.output);return
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1',NUMBA_CACHE_DIR='/private/tmp/m3_operator_pair_numba')
    summary={}
    for mode in ('legacy','ice_only'):
        out=ROOT/'bed_audit'/('paired_'+mode+'.json')
        try:
            result=subprocess.run([sys.executable,'-m','bed_audit.paired_development','--worker',mode,'--output',str(out)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=120)
            if result.returncode:
                summary[mode]=dict(status='failed',stderr=result.stderr[-3000:])
            else:
                data=json.loads(out.read_text());summary[mode]={k:data[k] for k in ('wall_s','loss_day','declared_minus_actual_melt_m3')};summary[mode]['status']='complete'
        except subprocess.TimeoutExpired:summary[mode]=dict(status='timeout',limit_s=120)
        print(mode,summary[mode],flush=True)
    (ROOT/'bed_audit'/'paired_summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':main()
