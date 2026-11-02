"""Safe selected text extraction from the byte-pinned public Kennicott archive."""
import argparse,hashlib,json,zipfile
from pathlib import Path,PurePosixPath
EXPECTED_BYTES=9948630
EXPECTED_MD5='f740fc1d0e6ff19c1b8049b427db0c3a'
EXPECTED_SHA256='a243a485eb83c51f54bf135cfecca9ca73e59ecc052fbe0ef94fbb70cca8044f'
PREFIX='Petersen_etal_2025_Kennicott_Glacier_Data/'
def extract(path,out):
 data=path.read_bytes()
 if len(data)!=EXPECTED_BYTES or hashlib.md5(data).hexdigest()!=EXPECTED_MD5 or hashlib.sha256(data).hexdigest()!=EXPECTED_SHA256:raise ValueError('unmatched original archive')
 records=[]
 with zipfile.ZipFile(path) as z:
  for member in z.infolist():
   if not member.filename.startswith(PREFIX):continue
   rel=PurePosixPath(member.filename[len(PREFIX):])
   if not rel.parts or rel.parts[0] not in ['Debris_Temperature_Profiles','Point_Mass_Balance','Debris_Thickness']:continue
   if rel.suffix.lower() not in ['.csv','.txt'] or any(p.startswith(('.', '~$')) or p=='..' for p in rel.parts):continue
   if member.is_dir() or member.file_size>5000000:continue
   raw=z.read(member);p=out.joinpath(*rel.parts);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
   records.append(dict(member=member.filename,path=str(rel),bytes=len(raw),zip_crc32=f'{member.CRC:08x}',sha256=hashlib.sha256(raw).hexdigest()))
 manifest=dict(source_record='Zenodo14625691 original verified archive, not assumed latest release',archive_bytes=len(data),archive_md5=EXPECTED_MD5,archive_sha256=hashlib.sha256(data).hexdigest(),members=records)
 (out.parent/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n');return manifest
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();m=extract(a.archive,a.out);print(len(m['members']),'verified members')
