"""Download only the official source archive, with hash/size checks."""
from pathlib import Path
import json,urllib.request,hashlib,os
R=Path(__file__).resolve().parent/'followup'

def main():
    meta=json.loads((R/'11581073_metadata.json').read_text());source=meta['files'][0]
    assert meta['metadata']['license']['id']=='cc-by-4.0'
    target=R/source['key'];partial=target.with_suffix('.download')
    if not target.exists():
        try:
            with urllib.request.urlopen(source['links']['self'],timeout=120) as response,partial.open('wb') as out:
                while True:
                    block=response.read(1024*1024)
                    if not block:break
                    out.write(block)
            data=partial.read_bytes()
            if len(data)!=source['size'] or 'md5:'+hashlib.md5(data).hexdigest()!=source['checksum']:raise ValueError('Source checksum/size mismatch')
            os.replace(partial,target)
        finally:
            if partial.exists():partial.unlink()
    data=target.read_bytes()
    assert len(data)==source['size'] and 'md5:'+hashlib.md5(data).hexdigest()==source['checksum']
    print('Official source verified:',target.name,len(data),'bytes')

if __name__=='__main__':main()
