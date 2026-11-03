"""Component checks; no outcome access or evaluation campaign."""
import json,os,re,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parent
runs=[]
for name in ['algorithm','alternative','alternative/neutral_flux','geometric','adaptive','novelty','turbulence']:
 p=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(R/name),'-p','test_*.py','-v'],capture_output=True,text=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1'),timeout=180)
 m=re.search(r'Ran (\d+) tests?',p.stderr);n=int(m.group(1))if m else 0
 runs.append(dict(component=name,tests=n,passed=p.returncode==0 and n>0,stdout=p.stdout,stderr=p.stderr));print(name,n,p.returncode,flush=True)
result=dict(tests=sum(x['tests']for x in runs),passed=all(x['passed']for x in runs),runs=runs,scope='Component verification, not independent field validation or historical novelty certification')
(R/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
if not result['passed']:raise SystemExit(1)
