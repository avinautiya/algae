"""Download the unchanged public source to a private directory; execute no PDF code."""
import argparse,hashlib
from pathlib import Path
from urllib.request import urlopen
import fitz
URL='https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C183D075F0E9ABCD7CE5251D6096FCA4/S0022143000206357a.pdf/outburst-flood-hazard-for-glacier-dammed-lac-de-rochemelon-france.pdf'
SHA='3ebd438caa17a6b6a8daf936cf228d528ac1fc65f73cffb529d0da68f86a507b'
p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True);a=p.parse_args();a.directory.mkdir(parents=True,exist_ok=True)
with urlopen(URL,timeout=60) as response:data=response.read(2000000)
if hashlib.sha256(data).hexdigest()!=SHA:raise ValueError('PDF changed; review source version before digitization')
path=a.directory/'vincent2010_publisher.pdf';path.write_bytes(data)
with fitz.open(path) as doc:
 d=doc.extract_image(91)
 if (d['width'],d['height'])!=(799,708):raise ValueError('Unexpected source image')
 (a.directory/'figure6_source.jpeg').write_bytes(d['image'])
print('Source verified; keep copyrighted source artifacts private unless rights permit distribution')
