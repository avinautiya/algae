"""Run bounded independent checks and authenticate preserved source blobs."""
import hashlib,json,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
EXPECTED={'validity_contract.py':'c63eeefa4b1841d4d7624138c054b19629017a1e','MODEL_VALIDITY_CONTRACT.md':'e84e83aba6d4b1e17328c10a8d9a04fa44171798','swiss_package_metadata_status.json':'67e29bd1c5572c2711fbbe78ceb54efd6ed65853','swiss_package_contract_result.json':'d42e54b8a6ad85cb127b435009e70a9b4fc004fb'}
def main():
 verified={}
 for name,expected in EXPECTED.items():
  data=(ROOT/'legacy'/name).read_bytes();actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
  verified[name]=dict(expected_git_blob=expected,actual_git_blob=actual,verified=expected==actual)
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 p=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT),'-p','test*.py','-v'],cwd=ROOT,env=env,capture_output=True,text=True)
 match=re.search(r'Ran (\d+) tests?',p.stderr)
 result=dict(passed=p.returncode==0 and all(v['verified'] for v in verified.values()),tests=int(match.group(1)) if match else None,pinned_repository_commit='05e7312e3bbf248acf9f646267a5b8ea06a8b32e',source_blobs=verified,stdout=p.stdout,stderr=p.stderr,scope='Typed metadata/migration and refusal behavior only; not environmental model validation')
 (ROOT/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['stdout','stderr']},indent=2))
 if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
