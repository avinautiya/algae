"""Run isolated component suites; no network, CFD or active-repo mutation."""
import json
from pathlib import Path
import re
import subprocess
import sys

def verify():
    root=Path(__file__).resolve().parent;results=[]
    for folder in ('data','geometry','forcing','transmission','integration'):
        suite=root/folder
        if not suite.is_dir():raise RuntimeError('missing component: '+folder)
        proc=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(suite),
                             '-p','test_*.py','-v'],text=True,capture_output=True,timeout=120)
        match=re.search(r'Ran (\d+) tests?',proc.stderr)
        count=int(match.group(1)) if match else 0
        results.append(dict(component=folder,tests=count,passed=proc.returncode==0 and count>0,
                            stdout=proc.stdout,stderr=proc.stderr))
        print(folder+': '+str(count)+' tests; '+('PASS' if results[-1]['passed'] else 'FAIL'),flush=True)
    result=dict(scope='component and conditional geometry checks; not environmental model validation',
                tests=sum(r['tests'] for r in results),passed=all(r['passed'] for r in results),results=results)
    (root/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['passed']:raise SystemExit(1)
    return result

if __name__=='__main__':verify()
