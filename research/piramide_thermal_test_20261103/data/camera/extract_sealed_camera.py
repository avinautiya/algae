"""Byte extraction/CRC only: does not deserialize any camera target."""
import json,struct,zlib,hashlib
from pathlib import Path
p=Path(__file__).parent
index=json.load(open('/private/tmp/melt-validation-20261103/model/environment/PIRAMIDE_ARCHIVE_INDEX.json'))
b=(p/'source_camera_range.bin').read_bytes();start=314012;members=[]
for r in index:
 if 'Ablation_Tower' not in r['name']:continue
 pos=r['offset']-start;h=struct.unpack_from('<4s5H3L2H',b,pos);q=pos+30+h[-2]+h[-1]
 raw=zlib.decompress(b[q:q+r['compressed_size']],-15)
 assert len(raw)==r['size'] and zlib.crc32(raw)==r['crc32']
 name=Path(r['name']).name;(p/name).write_bytes(raw)
 members.append({'file':name,'sha256':hashlib.sha256(raw).hexdigest(),'zip_crc32':r['crc32'],'size':len(raw),'outcomes_examined':False})
(p/'SEALED_SOURCE_MANIFEST.json').write_text(json.dumps({'source_doi':'10.5281/zenodo.11581073','members':members},indent=2)+'\n')
