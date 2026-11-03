"""Run bounded review regressions; never runs a recovery production campaign."""
import json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
 results=[]
 for name in ['endpoints','proof','scoring','benchmark','geometry']:
  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',NUMBA_CACHE_DIR=str(ROOT/'_runtime/numba_cache'))
  p=subprocess.run([sys.executable,'-m','unittest','discover','-p','test_*.py','-v'],cwd=ROOT/name,env=env,capture_output=True,text=True,timeout=240)
  m=re.search(r'Ran (\d+) tests?',p.stderr);n=int(m.group(1)) if m else 0
  results.append(dict(component=name,passed=p.returncode==0 and n>0,tests=n,stdout=p.stdout,stderr=p.stderr));print(name,results[-1]['passed'],n,flush=True)
 report=dict(passed=all(x['passed'] for x in results),tests=sum(x['tests'] for x in results),components=results,interpretation='bounded regression checks, not proof of all-error absence, field superiority or award outcomes')
 (ROOT/'VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
 if not report['passed']:raise SystemExit(1)
if __name__=='__main__':main()
