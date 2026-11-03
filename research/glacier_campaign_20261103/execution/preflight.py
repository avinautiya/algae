"""Bounded preflight on a verified COPY; never runs a glacier campaign."""
import argparse,hashlib,json,os,re,shutil,signal,subprocess,sys,time
from pathlib import Path
PIN='f1f85ee496dd39617fedf0f58202b943b3af3168'
PACKET=Path('research/cliff_repair_2026_11_01')

def verify_hashes(root):
 m=json.loads((root/'DELIVERY_MANIFEST.json').read_text());bad=[]
 for name,want in m['sha256'].items():
  p=Path(name)
  if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe manifest path')
  target=root/p
  if not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest()!=want:bad.append(name)
 if bad:raise ValueError('source/data manifest mismatch: '+', '.join(bad))
 return m

def check_scope(mode,covered):
 if mode not in ('physical_bare','conditional_quasi'):raise ValueError('unknown execution mode')
 if mode=='physical_bare' and covered:raise ValueError('physical covered-source mode is unsupported: drainage/storage closure required')
 return dict(mode=mode,covered=covered,environmental_validation=False,physical_covered_supported=False)

def bounded(cmd,cwd,seconds,log):
 if isinstance(seconds,bool) or seconds<=0:raise ValueError('positive timeout required')
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
 started=time.monotonic()
 with Path(log).open('x') as f:
  proc=subprocess.Popen(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
  try:code=proc.wait(timeout=seconds);timed_out=False
  except subprocess.TimeoutExpired:
   os.killpg(proc.pid,signal.SIGKILL);proc.wait();code=proc.returncode;timed_out=True
 return dict(command=cmd,returncode=code,timed_out=timed_out,wall_s=time.monotonic()-started,log=Path(log).name,threads_per_process=1)

GUARD_CODE='''import numpy as np
from phase.swe2d_v5 import Sim
for mode in ('solid','porous'):
 s=Sim(dx=.1,Lpool=.2,Lr=.4,Ltail=.2,Wl=.2,Q=.02,Tw=.5,ramp_ratio=0.)
 H=np.zeros((s.nx,s.ny));H[s.i0-1:s.i1,:s.nw]=.2
 s.set_deposit(H,mode);h=s.h.copy();z=s.zb.copy();ht=s.hT.copy()
 try:s.step_physical(.01)
 except ValueError as e:assert 'closure' in str(e) or 'covered' in str(e)
 else:raise AssertionError('unsupported covered physical source accepted')
 assert np.array_equal(h,s.h) and np.array_equal(z,s.zb) and np.array_equal(ht,s.hT)
print('both covered-source guards refuse before mutation')
'''

def run(repo,out,mode,covered,timeout,min_free_mb):
 repo=Path(repo).resolve();out=Path(out).resolve();check_scope(mode,covered)
 if out.exists():raise ValueError('output already exists; use a NEW directory')
 if repo==out or repo in out.parents:raise ValueError('output must be outside the checkout to protect frozen files')
 if not (repo/PACKET).is_dir():raise ValueError('repair packet missing; fetch the pinned repair branch')
 source=repo/PACKET;m=verify_hashes(source)
 for name in ['swe2d_v5.py','swe2d_v4_pinned.py']:
  installed=repo/'ice_cliff_survival/reactivation'/name
  if not installed.is_file() or installed.read_bytes()!=(source/'phase'/name).read_bytes():raise ValueError('installed/packet solver mismatch: '+name)
 free=shutil.disk_usage(out.parent if out.parent.exists() else repo).free
 if free<min_free_mb*1024**2:raise ValueError('insufficient free disk for verification; no cleanup attempted')
 out.mkdir(parents=True,exist_ok=False);work=out/'verified_copy'
 shutil.copytree(source,work,ignore=shutil.ignore_patterns('__pycache__','numba_cache','_runtime','*.pyc','*.nbi','*.nbc'))
 verify_hashes(work)
 git=subprocess.run(['git','-C',str(repo),'rev-parse','HEAD'],capture_output=True,text=True)
 result=dict(audited_source_pin=PIN,checkout_commit=git.stdout.strip() if git.returncode==0 else None,source_identity='packet hashes and installed solver copies verified; checkout identity is recorded separately',scope=check_scope(mode,covered),original_packet_unchanged=True,stages=[])
 try:
  result['stages'].append(bounded([sys.executable,'verify_repairs.py'],work,timeout,out/'tests.log'))
  if result['stages'][-1]['returncode']!=0:raise ValueError('repair tests failed or timed out')
  v=json.loads((work/'VERIFICATION.json').read_text())
  if v.get('passed') is not True or v.get('tests')!=66:raise ValueError('expected all66 pinned tests; source update needs a new reviewed freeze')
  result['tests_passed']=v['tests']
  result['stages'].append(bounded([sys.executable,'-c',GUARD_CODE],work,timeout,out/'covered_guard.log'))
  if result['stages'][-1]['returncode']!=0:raise ValueError('covered-source refusal check failed')
  verify_hashes(source)
  result['preflight_passed']=True
 except Exception as e:
  result.update(preflight_passed=False,error=str(e))
 (out/'PREFLIGHT_RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
 return result

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--mode',choices=['physical_bare','conditional_quasi'],required=True);p.add_argument('--covered',action='store_true');p.add_argument('--timeout-s',type=int,default=600);p.add_argument('--min-free-mb',type=int,default=256);a=p.parse_args()
 try:r=run(a.repo,a.out,a.mode,a.covered,a.timeout_s,a.min_free_mb)
 except (OSError,ValueError) as e:r=dict(preflight_passed=False,error=str(e))
 print(json.dumps(r,indent=2));raise SystemExit(0 if r['preflight_passed'] else 2)
if __name__=='__main__':main()
