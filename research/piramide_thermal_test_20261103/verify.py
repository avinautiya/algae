import json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
    results=[]
    for folder in ['data','data/camera','model','analysis','review']:
        p=subprocess.run([sys.executable,'-m','unittest','discover','-p','test_*.py','-v'],cwd=ROOT/folder,
            env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1'),capture_output=True,text=True,timeout=180)
        match=re.search(r'Ran (\d+) tests?',p.stderr);count=int(match.group(1)) if match else 0
        results.append(dict(component=folder,tests=count,passed=p.returncode==0 and count>0,stdout=p.stdout,stderr=p.stderr))
        print(folder,count,p.returncode,flush=True)
    report=dict(tests=sum(r['tests'] for r in results),passed=all(r['passed'] for r in results),results=results,
        scope='Loader, inference, numerical and sealing checks; no superiority certification')
    (ROOT/'VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['passed']:raise SystemExit(1)
if __name__=='__main__':main()
