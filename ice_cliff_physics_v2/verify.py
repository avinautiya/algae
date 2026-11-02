"""Run isolated mathematical/implementation checks; no field accuracy inference."""
import json,os,re,subprocess,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def main():
    groups=[('.', 'test_model.py'),('.', 'test_upgrades.py'),('energy_ports','test_*.py'),
            ('thermal','test_*.py'),('transfer','test_*.py'),('face_geometry','test_*.py'),('mechanism','test_*.py')]
    results=[]
    for folder,pattern in groups:
        env=os.environ.copy();env['PYTHONPATH']=os.pathsep.join([str(ROOT/folder),str(ROOT)])
        process=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT/folder),'-t',str(ROOT),'-p',pattern],
            capture_output=True,text=True,env=env,timeout=120)
        output=process.stdout+process.stderr;match=re.search(r'Ran (\d+) tests?',output)
        results.append(dict(group=folder,pattern=pattern,exit_code=process.returncode,
                            tests=int(match.group(1)) if match else 0,output=output))
        print(folder,pattern,'PASS' if process.returncode==0 else 'FAIL',flush=True)
    result=dict(scope='analytic, conservation, numerical and implementation checks; not field validation',
                passed=all(r['exit_code']==0 for r in results),tests=sum(r['tests'] for r in results),groups=results)
    (ROOT/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit(1)
    print(json.dumps(dict(passed=True,tests=result['tests'])),flush=True)

if __name__=='__main__':main()
