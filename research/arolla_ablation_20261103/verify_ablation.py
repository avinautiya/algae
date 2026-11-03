import json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
 results=[]
 for name in ['targets','reproduction','model','timing','model/environment','reproduction/robustness']:
  p=subprocess.run([sys.executable,'-m','unittest','discover','-p','test_*.py','-v'],cwd=ROOT/name,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1'),capture_output=True,text=True,timeout=180)
  m=re.search(r'Ran (\d+) tests?',p.stderr);n=int(m.group(1)) if m else 0
  results.append(dict(component=name,tests=n,passed=p.returncode==0 and n>0,stdout=p.stdout,stderr=p.stderr));print(name,n,p.returncode,flush=True)
 r=dict(passed=all(x['passed'] for x in results),tests=sum(x['tests'] for x in results),results=results,scope='unit, loader and bounded-arithmetic checks, not ablation-performance certification')
 (ROOT/'VERIFICATION.json').write_text(json.dumps(r,indent=2)+'\n')
 if not r['passed']:raise SystemExit(1)
if __name__=='__main__':main()
