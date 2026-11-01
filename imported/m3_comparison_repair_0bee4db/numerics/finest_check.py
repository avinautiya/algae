"""Plan or isolate one finest-grid check; never launches without --execute.

The pinned source is EXPECTED TO REFUSE the geometry/latent preflight. This
runner must not be used to manufacture a finest-grid endorsement of bad code.
"""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import types
import numpy as np

PIN='f3602e35e1595bcdbee4208a354a34944f1ae974'


def load_snapshot(path,cache_dir):
    os.environ['NUMBA_CACHE_DIR']=str(cache_dir)
    cache_dir.mkdir(parents=True,exist_ok=True)
    path=path.resolve()
    spec=json.loads((path/'frozen/production_spec.json').read_text())
    actual={name:hashlib.sha256((path/name).read_bytes()).hexdigest() for name in ('swe2d.py','cases.py')}
    if any(actual[k]!=spec['solver']['file_sha256'][k] for k in actual):
        raise ValueError('snapshot differs from its frozen solver hashes')
    package=types.ModuleType('_m3_finest_snapshot');package.__path__=[str(path)]
    sys.modules[package.__name__]=package
    return importlib.import_module(package.__name__+'.cases'),actual


def operator_preflight(cases):
    """Independent actual-geometry vs reported-volume checks; no hydrodynamics."""
    s=cases.build('strip',dx=.1)
    s.last_qacc=np.zeros_like(s.zb)
    before=s.zb.copy();reported=s.lower_bed(0.)
    changed=np.abs(s.zb-before)>1e-14
    zero=dict(nonice_changed_cells=int(np.count_nonzero(changed&~s.ice)),
              total_changed_cells=int(np.count_nonzero(changed)),
              maximum_upward_change_m=float(np.max(s.zb-before)),
              reported_volume_m3=float(reported[0]))
    # Separate fresh object: force one ice column beyond its available thickness.
    s=cases.build('strip',dx=.1);s.last_qacc=np.zeros_like(s.zb)
    i,j=np.argwhere(s.ice)[0]
    available=max(float(s.zb[i,j]-s.pool_floor),0.)
    proposed=available+.1
    s.last_qacc[i,j]=proposed*cases.RHO_I*cases.LF
    before=s.zb.copy();reported=s.lower_bed(1.)
    actual=float(((before-s.zb)*s.ice).sum()*s.dx**2)
    cap=dict(reported_volume_m3=float(reported[0]),actual_ice_lowering_m3=actual,
             reported_minus_actual_m3=float(reported[0]-actual),
             forced_proposed_lowering_m=proposed)
    passed=zero['total_changed_cells']==0 and abs(zero['reported_volume_m3'])<1e-14 and abs(cap['reported_minus_actual_m3'])<1e-12
    return dict(passed=passed,zero_dt=zero,cap_activation=cap,
                scope='geometry/volume operator checks only; not the full verification gate')


def plan(days,removal):
    # E04 measured medium21d2028s; using Ncells ~h^-2 and CFL steps ~h^-1.
    est=2028.*8*days/21.
    return dict(dx_m=.0125,days=days,removal_day=removal,workers=1,
                medium_case_reference_wall_s=2028.,
                estimated_compute_wall_s=est,
                estimate_scope='single-worker cubic grid-cost extrapolation; excludes uncertain startup, contention and changed toolchain',
                supports=('early contact-loss refinement only' if removal is None else 'declared removal case over the actually simulated horizon'),
                does_not_support='short early-loss run does not establish day12 recovery or21d melt convergence',
                checkpoint_support=False,
                checkpoint_note='stock cases.run has no checkpoint API; timeout preserves progress metadata, not a resumable state')


def atomic_json(path,data):
    temp=path.with_name(path.name+'.tmp');temp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');temp.replace(path)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot',type=Path,default=Path(__file__).resolve().parents[1]/'snapshot')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--days',type=float,default=3.5)
    parser.add_argument('--removal-day',type=float)
    parser.add_argument('--wall-budget-s',type=float,default=7200.)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    if not np.isfinite(args.days) or args.days<=0 or not np.isfinite(args.wall_budget_s) or args.wall_budget_s<=0:parser.error('positive finite horizon/budget required')
    if args.removal_day is not None and (not np.isfinite(args.removal_day) or not 0<=args.removal_day<=args.days):parser.error('removal outside horizon')
    if not np.isclose(args.days*4,round(args.days*4),rtol=0,atol=1e-10):parser.error('stock runner requires horizon aligned to6h')
    if args.removal_day is not None and not np.isclose(args.removal_day*4,round(args.removal_day*4),rtol=0,atol=1e-10):parser.error('stock runner requires removal aligned to6h')
    # Refuse to mix snapshot source and outputs or overwrite an existing release.
    args.output=args.output.resolve();snapshot=args.snapshot.resolve()
    if args.output==snapshot or snapshot in args.output.parents:parser.error('output must be outside immutable snapshot')
    args.output.mkdir(parents=True,exist_ok=True)
    if not args.worker and (args.output/'result.json').exists():parser.error('existing result will not be overwritten')
    cases,hashes=load_snapshot(snapshot,args.output/'numba_cache')
    gate=operator_preflight(cases)
    document=dict(reference_report_commit=PIN,source_hashes=hashes,
                  source_identification='content hashes; report pin does not authenticate arbitrary --snapshot directories',
                  plan=plan(args.days,args.removal_day),preflight=gate,
                  status='planned' if gate['passed'] else 'blocked_geometry_operator')
    atomic_json(args.output/'plan.json',document)
    if not args.execute:return 0
    if not gate['passed']:
        print('REFUSED: geometry/latent preflight failed; repair and refreeze before finest run',file=sys.stderr);return 2
    if args.worker:
        atomic_json(args.output/'progress.json',dict(status='running',simulated_time_progress_available=False,checkpoint_support=False))
        started=time.monotonic()
        result=cases.run('strip',dx=.0125,days=args.days,removal_day=args.removal_day,dtm=21600.)
        atomic_json(args.output/'result.json',result)
        atomic_json(args.output/'progress.json',dict(status='complete',wall_s=time.monotonic()-started,simulated_days=args.days))
        return 0
    command=[sys.executable,str(Path(__file__).resolve()),'--snapshot',str(snapshot),'--output',str(args.output),'--days',str(args.days),'--execute','--worker']
    if args.removal_day is not None:command+=['--removal-day',str(args.removal_day)]
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):env[key]='1'
    with (args.output/'worker.log').open('w') as log:
        proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        try:code=proc.wait(timeout=args.wall_budget_s)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM)
            try:proc.wait(timeout=10.)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
            atomic_json(args.output/'progress.json',dict(status='wall_budget_terminated',budget_s=args.wall_budget_s,checkpoint_support=False))
            return 3
    return code


if __name__=='__main__':sys.exit(main())
