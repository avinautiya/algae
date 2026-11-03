import json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
    results=[]
    for name in ['theory','physics','literature','moisture','moisture/review','data/moisture']:
        p=subprocess.run([sys.executable,'-m','unittest','discover','-p','test_*.py','-v'],cwd=ROOT/name,
            capture_output=True,text=True,timeout=180,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1'))
        match=re.search(r'Ran (\d+) tests?',p.stderr);n=int(match.group(1))if match else 0
        results.append(dict(component=name,tests=n,passed=p.returncode==0 and n>0,stdout=p.stdout,stderr=p.stderr));print(name,n,p.returncode,flush=True)
    r=dict(tests=sum(x['tests']for x in results),passed=all(x['passed']for x in results),results=results,
        scope='Mathematical/numerical/source checks; not new field validation or historical originality certification')
    (ROOT/'VERIFICATION.json').write_text(json.dumps(r,indent=2)+'\n')
    if not r['passed']:raise SystemExit(1)
if __name__=='__main__':main()
