"""Reproduce isolated repair tests; never launches a glacier campaign."""
import json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
 rows=[]
 for group in ['phase','contract','evidence','roi']:
  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',NUMBA_CACHE_DIR=str(ROOT/'_runtime/numba_cache'),PYTHONPATH=str(ROOT/group))
  p=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/group),'-p','test_*.py','-v'],cwd=ROOT/group,env=env,text=True,capture_output=True)
  m=re.search(r'Ran (\d+) tests?',p.stderr)
  rows.append(dict(group=group,passed=p.returncode==0,tests=int(m.group(1)) if m else None,stdout=p.stdout,stderr=p.stderr))
  print(group,'PASS' if p.returncode==0 else 'FAIL',rows[-1]['tests'],flush=True)
 result=dict(passed=all(r['passed'] for r in rows),tests=sum(r['tests'] or 0 for r in rows),groups=rows,interpretation='software/source-step/numerical checks only; not field validation or novelty certification')
 (ROOT/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
 if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
