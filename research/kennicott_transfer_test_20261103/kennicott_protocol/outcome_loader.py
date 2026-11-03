"""Role-gated native net-ablation loader. No import-time network or outcome reads.
Publication/role receipts are an auditable logical gate, not cryptographic secrecy:
TRAIN and EVALUATION rows share a public CSV. Only authorized rows are converted.
"""
import csv,datetime as dt,hashlib,io,json,math,struct,urllib.request,zlib
from pathlib import Path
SOURCE_SHA256='4b1fd8e9d7d0a9aa7365bfc40d3e90ac70948361aaa4cccb993412f2b39e8cd8'
SOURCE_CRC32=1387455814
SOURCE_BYTES=30661
ARCHIVE_URL='https://zenodo.org/records/16389806/files/Petersen_etal_2025_Kennicott_Glacier_Data.zip?download=1'
MEMBER='Petersen_etal_2025_Kennicott_Glacier_Data/Point_Mass_Balance/point_mass_balance.csv'
MEMBER_OFFSET=4216008
MEMBER_COMPRESSED_BYTES=6254
OPERATOR='negative_source_mass_change'
BINDINGS={'train':[{'site':'P_AWS_1','start':'2020-07-25','end':'2020-09-06','period_days':43,'debris_mm':210}], 'evaluation':[{'site':s,'start':'2023-06-08','end':'2023-09-02','period_days':86,'debris_mm':d} for s,d in [('P11',240),('P12',290),('P13',255)]]}

def release_plan():
 return {'approved_operator':OPERATOR,'unit':'mmwe','source_csv_sha256':SOURCE_SHA256,'source_member_crc32':SOURCE_CRC32,'source_member':MEMBER,'bindings':json.loads(json.dumps(BINDINGS)),'target_interpretation':'conditional_net_ice_ablation_loss_not_verified_gross_melt','logical_public_source_gate_not_cryptographic_blind':True}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def _path(rpath,p):
 q=Path(p);return q if q.is_absolute() else Path(rpath).parent/q

def verify_receipt(receipt_path,*,expected_role):
 if expected_role not in BINDINGS:raise ValueError('Explicit authorized role required')
 rp=Path(receipt_path);r=json.loads(rp.read_text())
 if r.get('schema')!='kennicott-outcome-release-v1' or r.get('access_role')!=expected_role:raise ValueError('Wrong receipt schema or access role')
 if r.get('freeze_published') is not True or r.get('predictions_published') is not True:raise ValueError('Published spec and predictions required')
 if not isinstance(r.get('publication_reference'),str) or not r['publication_reference'].strip():raise ValueError('Publication reference required')
 if r.get('prediction_stage')!={'train':'normalized_training','evaluation':'final_evaluation'}[expected_role]:raise ValueError('Wrong prediction release stage')
 for pkey,hkey in [('frozen_spec_path','frozen_spec_sha256'),('predictions_path','predictions_sha256')]:
  if not isinstance(r.get(pkey),str) or not r[pkey] or not isinstance(r.get(hkey),str) or sha(_path(rp,r[pkey]))!=r[hkey]:raise ValueError('Frozen artifact hash mismatch')
 spec=json.loads(_path(rp,r['frozen_spec_path']).read_text())
 if spec.get('native_outcome_release')!=release_plan():raise ValueError('Native source, operator and roles not bound into frozen spec')
 if r.get('approved_operator')!=OPERATOR:raise ValueError('Unapproved operator')
 return r

def _range(start,end):
 req=urllib.request.Request(ARCHIVE_URL,headers={'Range':f'bytes={start}-{end}'})
 with urllib.request.urlopen(req,timeout=45) as r:
  if r.status!=206:raise ValueError('Server refused bounded range; full archive not read')
  b=r.read(end-start+1)
 if len(b)!=end-start+1:raise ValueError('Truncated source byte range')
 return b

def _fetch_verified_native():
 h=_range(MEMBER_OFFSET,MEMBER_OFFSET+29);f=struct.unpack('<4s5H3L2H',h)
 if f[0]!=b'PK\x03\x04' or f[3]!=8:raise ValueError('Unexpected ZIP local-header compression')
 name_extra=_range(MEMBER_OFFSET+30,MEMBER_OFFSET+30+f[-2]+f[-1]-1)
 if name_extra[:f[-2]].decode('utf-8')!=MEMBER:raise ValueError('Unexpected ZIP member')
 start=MEMBER_OFFSET+30+f[-2]+f[-1];compressed=_range(start,start+MEMBER_COMPRESSED_BYTES-1);raw=zlib.decompress(compressed,-15)
 if len(raw)!=SOURCE_BYTES or zlib.crc32(raw)&0xffffffff!=SOURCE_CRC32 or hashlib.sha256(raw).hexdigest()!=SOURCE_SHA256:raise ValueError('Native outcome source hash/CRC mismatch')
 return raw

def _date(value):
 # Source CSV uses US month/day/year despite metadata describing ISO dates.
 return dt.datetime.strptime(value,'%m/%d/%Y').date().isoformat()

def _select_rows(raw,bindings):
 rows=list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))));required={'site','date_start','date_end','period_length','mass_change_kg_m2','debris_thick_mm','surface_type_start','surface_type_end','seasonal_snow_depth_start_mm','seasonal_snow_depth_end_mm','summer_snow_depth_end_mm'}
 if not rows or not required.issubset(rows[0]):raise ValueError('Native stake schema changed')
 output=[]
 for binding in bindings:
  matches=[]
  for row in rows:
   if row['site']!=binding['site']:continue
   try:key=(_date(row['date_start']),_date(row['date_end']))
   except ValueError:continue
   if key==(binding['start'],binding['end']):matches.append(row)
  if len(matches)!=1:raise ValueError('Authorized native row absent or duplicated')
  row=matches[0]
  if int(row['period_length'])!=binding['period_days'] or (dt.date.fromisoformat(binding['end'])-dt.date.fromisoformat(binding['start'])).days!=binding['period_days']:raise ValueError('Native period mismatch')
  if float(row['debris_thick_mm'])!=binding['debris_mm']:raise ValueError('Stake geometry changed')
  if any(row[k].strip().lower()!='debris' for k in ['surface_type_start','surface_type_end']):raise ValueError('Unsupported native surface component')
  if any(float(row[k])!=0 for k in ['seasonal_snow_depth_start_mm','seasonal_snow_depth_end_mm']):raise ValueError('Native snow cover violates frozen target interpretation')
  if row['summer_snow_depth_end_mm'].strip() and float(row['summer_snow_depth_end_mm'])!=0:raise ValueError('Native summer snow violates frozen target interpretation')
  # First numerical mass conversion occurs here, and only for selected role rows.
  source=float(row['mass_change_kg_m2'])
  if not math.isfinite(source):raise ValueError('Nonfinite authorized native target')
  target=-source
  output.append({**binding,'observed_net_ablation_mmwe':target,'negative_net_ablation_retained':target<0,'operator':OPERATOR,'density_multiplier_applied':False,'native_unit':'kg_m2_equivalent_mmwe','native_source_field':'mass_change_kg_m2','conditional_interpretation':'net loss; not independently verified gross phase-change melt','uncertainty':'No statistical measurement interval inferred from source rounding or sample count'})
 return output

def load_after_freeze(receipt_path,*,expected_role):
 receipt=verify_receipt(receipt_path,expected_role=expected_role)
 raw=_fetch_verified_native() # Never called before both publication/role checks.
 return {'access_role':expected_role,'receipt_sha256':sha(receipt_path),'spec_sha256':receipt['frozen_spec_sha256'],'predictions_sha256':receipt['predictions_sha256'],'source_csv_sha256':SOURCE_SHA256,'source_crc32':SOURCE_CRC32,'rows':_select_rows(raw,BINDINGS[expected_role]),'release_scope':'Only exact native rows for authorized role converted and returned; no complete outcome table exported.'}
