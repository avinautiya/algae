"""Run independent field-component checks without re-fitting on held-out data."""
import json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 groups=['solver','protocol','ablation_benchmark','temperature_benchmark','layered_solver']
 records=[]
 for name in groups:
  folder=ROOT/'field_thermal'/name
  if not any(folder.glob('test*.py')):
   records.append(dict(group=name,status='not_present'));continue
  env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',PYTHONPATH=str(folder))
  proc=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(folder),'-p','test*.py','-v'],cwd=ROOT,env=env,text=True,capture_output=True)
  records.append(dict(group=name,status='passed' if proc.returncode==0 else 'failed',stdout=proc.stdout,stderr=proc.stderr))
  print(name,records[-1]['status'],flush=True)
 passed=all(x['status']=='passed' for x in records)
 (ROOT/'FIELD_VERIFICATION.json').write_text(json.dumps(dict(passed=passed,groups=records),indent=2)+'\n')
 if not passed:raise SystemExit(1)
if __name__=='__main__':main()
