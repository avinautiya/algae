"""Check every file in the published evidence manifest."""
from pathlib import Path
import hashlib,json
R=Path(__file__).resolve().parent

def main():
    manifest=json.loads((R/'MANIFEST.json').read_text());errors=[]
    for row in manifest['files']:
        p=R/row['path']
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:errors.append(row['path'])
    if errors:raise SystemExit('Missing or modified artifacts: '+', '.join(errors))
    print('PASS:',len(manifest['files']),'published file hashes verified')

if __name__=='__main__':main()
