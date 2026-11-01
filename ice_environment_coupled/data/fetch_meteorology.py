"""Retrieve a pinned published Lirung AWS example, preserving original bytes."""
import hashlib,json,urllib.request
from pathlib import Path
COMMIT='ea5ca8bae9b3c67ec3cf2a8b56c724d87bb1909c'
BLOB='db18c3c75a55ed1642f2c91c9aac2001afc10a7c'
URL='https://raw.githubusercontent.com/pburi/CliffEBM/'+COMMIT+'/Meteodata/meteodata_AWS2013.txt'
def fetch(root):
    with urllib.request.urlopen(URL,timeout=30) as response:data=response.read(1_000_001)
    if len(data)>1_000_000:raise RuntimeError('unexpected source size')
    actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if actual!=BLOB:raise RuntimeError('source Git blob mismatch')
    root.mkdir(parents=True,exist_ok=True);(root/'meteodata_AWS2013.txt').write_bytes(data)
    manifest=dict(source_url=URL,source_commit=COMMIT,source_git_blob=BLOB,git_blob_verified=True,
                  sha256=hashlib.sha256(data).hexdigest(),bytes=len(data),
                  attribution='Pascal Buri, CliffEBM Lirung example. Public repository GPL-3.0; data-specific reuse terms must be confirmed before redistribution.',
                  scope='Atmospheric station forcing, not stream runoff, water temperature or a measured cliff-face radiation field.')
    (root/'METEOROLOGY_PROVENANCE.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__':fetch(Path(__file__).parent)
