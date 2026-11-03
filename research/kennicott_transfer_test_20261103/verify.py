import json,os,re,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parent;r=[]
for name in ['ROOT','continuous_goal','initial_state','algorithm','review','kennicott_model','kennicott_forcing','kennicott_protocol']:
 cmd=[sys.executable,'-m','unittest',name+'.test_boundary'if name=='kennicott_forcing'else name+'.test_outcome_loader','-v']if name in ['kennicott_forcing','kennicott_protocol']else[sys.executable,'-m','unittest','discover','-s',str(R if name=='ROOT'else R/name),'-p','test_locked_experiment.py'if name=='ROOT'else'test_*.py','-v']
 p=subprocess.run(cmd,cwd=R,capture_output=True,text=True,timeout=180,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1'))
 m=re.search(r'Ran (\d+) tests?',p.stderr);n=int(m.group(1))if m else 0
 r.append(dict(component=name,tests=n,passed=p.returncode==0 and n>0,stderr=p.stderr,stdout=p.stdout));print(name,n,p.returncode,flush=True)
x=dict(passed=all(v['passed']for v in r),tests=sum(v['tests']for v in r),results=r,scope='Numerical/component/source-gate verification, not field accuracy or historical originality')
(R/'VERIFICATION.json').write_text(json.dumps(x,indent=2)+'\n')
if not x['passed']:raise SystemExit(1)
