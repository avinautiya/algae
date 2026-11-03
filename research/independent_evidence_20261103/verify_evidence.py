import json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
 out=[]
 for folder in ['reproduction','extension','decision','coverage']:
  p=subprocess.run([sys.executable,'-m','unittest','discover','-p','test_*.py','-v'],cwd=ROOT/folder,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1'),capture_output=True,text=True,timeout=180)
  m=re.search(r'Ran (\d+) tests?',p.stderr);n=int(m.group(1)) if m else 0
  out.append(dict(component=folder,tests=n,passed=p.returncode==0 and n>0,stdout=p.stdout,stderr=p.stderr));print(folder,n,p.returncode,flush=True)
 r=dict(passed=all(d['passed'] for d in out),tests=sum(d['tests'] for d in out),components=out,scope='unit and synthetic checks only; held-out signals not scored by verifier')
 (ROOT/'VERIFICATION.json').write_text(json.dumps(r,indent=2)+'\n')
 if not r['passed']:raise SystemExit(1)
if __name__=='__main__':main()
