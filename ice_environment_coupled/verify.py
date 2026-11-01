"""Bounded, network-free checks; no CFD or active source changes."""
import hashlib,json,re,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
commands=[['-m','unittest','-v','energy_ports.test_energy'],
          ['-m','unittest','discover','-s',str(root),'-p','test_model.py','-v'],
          ['-m','unittest','discover','-s',str(root/'data'),'-p','test_*.py','-v']]
def verify():
    hashes=json.loads((root/'vendor/PROVENANCE.json').read_text())
    for name,record in hashes.items():
        if hashlib.sha256((root/'vendor'/name).read_bytes()).hexdigest()!=record['sha256']:
            raise RuntimeError('vendored source changed: '+name)
    results=[]
    for args in commands:
        p=subprocess.run([sys.executable,*args],cwd=root,capture_output=True,text=True,timeout=120)
        match=re.search(r'Ran (\d+) tests?',p.stderr)
        count=int(match.group(1)) if match else 0
        results.append(dict(command=args,tests=count,passed=p.returncode==0 and count>0,stdout=p.stdout,stderr=p.stderr))
        print(str(count)+' tests: '+('PASS' if results[-1]['passed'] else 'FAIL'),flush=True)
    report=dict(tests=sum(r['tests'] for r in results),passed=all(r['passed'] for r in results),
                vendor_hashes_verified=True,scope='component/integration checks, not environmental validation',results=results)
    (root/'VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['passed']:raise SystemExit(1)
if __name__=='__main__':verify()
