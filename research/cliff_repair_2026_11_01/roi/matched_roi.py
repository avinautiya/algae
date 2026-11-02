"""Float64 geometric melt observation on an explicitly fixed rectangular ROI.

The output manifest verifies byte integrity, not who ran a simulation or whether
its inputs are environmental truth. No legacy output is silently authenticated.
"""
from dataclasses import dataclass
import hashlib,io,json,math,os,re,tempfile,zipfile
from pathlib import Path
import numpy as np

def digest_bytes(b):return hashlib.sha256(b).hexdigest()
def array_sha(a):
 a=np.ascontiguousarray(a);return digest_bytes(str(a.shape).encode()+a.dtype.str.encode()+a.tobytes())
def number(x):
 try:return type(x) in (int,float) and math.isfinite(x)
 except (OverflowError,ValueError):return False
def canonical(x):return json.dumps(x,sort_keys=True,allow_nan=False,separators=(',',':')).encode()

def strict_json(raw):
 def bad(x):raise ValueError('nonfinite JSON constant '+x)
 def unique(pairs):
  out={}
  for k,v in pairs:
   if k in out:raise ValueError('duplicate JSON key '+k)
   out[k]=v
  return out
 return json.loads(raw,parse_constant=bad,object_pairs_hook=unique)

def valid_sha(x):return isinstance(x,str) and re.fullmatch('[0-9a-f]{64}',x) is not None and len(set(x))>1

def validate_record(r,*,mask=None,initial=None,current=None):
 keys={'roi_id','bounds_m','coordinate_frame','area_m2','dx_m','physical_time_s','geometric_melt_m3','heat_implied_melt_m3','geometry_heat_residual_m3','absolute_ledger_tolerance_m3','mask_sha256','initial_bed_sha256','current_bed_sha256','provenance','scope'}
 if not isinstance(r,dict) or set(r)!=keys:raise ValueError('record schema mismatch')
 canonical(r)
 if not all(isinstance(r[k],str) and r[k].strip() for k in ['roi_id','coordinate_frame','scope']):raise ValueError('ROI/frame/scope required')
 if not isinstance(r['provenance'],dict) or not r['provenance']:raise ValueError('ROI provenance required')
 if not all(number(r[k]) and r[k]>=0 for k in ['physical_time_s','geometric_melt_m3','absolute_ledger_tolerance_m3']):raise ValueError('invalid time/volume/tolerance')
 if not all(number(r[k]) and r[k]>0 for k in ['dx_m','area_m2']):raise ValueError('positive finite grid/area required')
 if not all(valid_sha(r[k]) for k in ['mask_sha256','initial_bed_sha256','current_bed_sha256']):raise ValueError('invalid array hashes')
 bounds=r['bounds_m'];dx=r['dx_m']
 if not isinstance(bounds,list) or len(bounds)!=4 or not all(number(x) for x in bounds):raise ValueError('invalid ROI bounds')
 x0,x1,y0,y1=bounds
 if not(0<=x0<x1 and 0<=y0<y1) or any(abs(v/dx-round(v/dx))>1e-9 for v in bounds):raise ValueError('invalid/unaligned ROI bounds')
 if abs((x1-x0)*(y1-y0)-r['area_m2'])>1e-12*max(1,r['area_m2']):raise ValueError('bounds/area mismatch')
 heat=r['heat_implied_melt_m3'];res=r['geometry_heat_residual_m3']
 if heat is None:
  if res is not None:raise ValueError('residual without heat reference')
 else:
  if not number(heat) or heat<0 or not number(res):raise ValueError('invalid heat volume/residual')
  expected_residual=r['geometric_melt_m3']-heat
  if res!=expected_residual or abs(expected_residual)>r['absolute_ledger_tolerance_m3']:raise ValueError('heat/geometric ledger mismatch')
 if mask is not None:
  if mask.dtype!=np.bool_ or mask.ndim!=2 or min(mask.shape)<3 or not mask.any():raise ValueError('invalid mask')
  if x1>(mask.shape[0]-2)*dx+1e-12 or y1>(mask.shape[1]-2)*dx+1e-12:raise ValueError('bounds outside data grid')
  expected=np.zeros(mask.shape,dtype=bool);expected[1+round(x0/dx):1+round(x1/dx),1+round(y0/dx):1+round(y1/dx)]=True
  if not np.array_equal(mask,expected) or array_sha(mask)!=r['mask_sha256']:raise ValueError('bounds/mask mismatch')
  if any(x is None or x.dtype!=np.float64 or x.shape!=mask.shape or not np.isfinite(x).all() for x in [initial,current]):raise ValueError('matching finite float64 beds required')
  if array_sha(initial)!=r['initial_bed_sha256'] or array_sha(current)!=r['current_bed_sha256']:raise ValueError('bed hash mismatch')
  if (current[mask]-initial[mask]>1e-12).any():raise ValueError('unsupported ice growth')
  volume=float(np.sum((initial-current)[mask])*dx**2)
  if not math.isfinite(volume) or abs(volume-r['geometric_melt_m3'])>1e-12*max(1,abs(volume)):raise ValueError('geometric volume mismatch')
 return r

@dataclass(frozen=True)
class Reference:
 bed_m:np.ndarray
 mask:np.ndarray
 dx_m:float
 bounds_m:tuple
 frame:str
 roi_id:str
 provenance:dict

def capture(sim,*,roi_id,bounds_m,frame,provenance):
 """Bounds(x0,x1,y0,y1) use cell corners, origin0 for interior Sim cells.

 Refuse nonaligned bounds instead of silently moving the scored footprint.
 This helper supports only the current uniform Cartesian grid convention.
 """
 if not isinstance(roi_id,str) or not roi_id.strip() or not isinstance(frame,str) or not frame.strip():raise ValueError('explicit ROI and coordinate frame required')
 if not isinstance(provenance,dict) or not provenance:raise ValueError('source-labelled ROI provenance required')
 if not isinstance(bounds_m,(list,tuple)) or len(bounds_m)!=4 or not all(number(x) for x in bounds_m):raise ValueError('four finite bounds required')
 x0,x1,y0,y1=bounds_m;dx=float(sim.dx)
 if not number(dx) or dx<=0 or not(0<=x0<x1<=sim.nx*dx and 0<=y0<y1<=sim.ny*dx):raise ValueError('ROI outside grid')
 if any(abs(v/dx-round(v/dx))>1e-9 for v in bounds_m):raise ValueError('ROI boundaries must align with cell edges')
 shape=(sim.nx+2,sim.ny+2)
 z=np.asarray(sim.zb,dtype=np.float64)
 if z.shape!=shape or not np.isfinite(z).all():raise ValueError('finite padded bed grid required')
 mask=np.zeros(shape,dtype=bool);mask[1+round(x0/dx):1+round(x1/dx),1+round(y0/dx):1+round(y1/dx)]=True
 if not np.all(np.asarray(sim.ice,dtype=bool)[mask]):raise ValueError('requested rectangular footprint includes nonice cells')
 bed=z.copy();mask=mask.copy();bed.flags.writeable=False;mask.flags.writeable=False
 return Reference(bed,mask,dx,tuple(bounds_m),frame,roi_id,json.loads(canonical(provenance)))

def observe(sim,ref,*,physical_time_s,heat_implied_volume_m3=None,absolute_ledger_tolerance_m3=1e-10):
 if not number(physical_time_s) or physical_time_s<0:raise ValueError('finite nonnegative physical time required')
 if not number(absolute_ledger_tolerance_m3) or absolute_ledger_tolerance_m3<0:raise ValueError('finite nonnegative tolerance required')
 if float(sim.dx)!=ref.dx_m:raise ValueError('grid changed during observation')
 bed=np.asarray(sim.zb,dtype=np.float64)
 if bed.shape!=ref.bed_m.shape or not np.isfinite(bed).all():raise ValueError('bed shape/values invalid')
 dz=(ref.bed_m-bed)[ref.mask]
 if (dz < -1e-12).any():raise ValueError('ice growth/refreezing outside supported observation')
 volume=float(np.sum(dz)*ref.dx_m**2)
 if not math.isfinite(volume):raise ValueError('volume overflows')
 residual=None
 if heat_implied_volume_m3 is not None:
  if not number(heat_implied_volume_m3) or heat_implied_volume_m3<0:raise ValueError('finite nonnegative heat-implied volume required')
  residual=volume-heat_implied_volume_m3
  if abs(residual)>absolute_ledger_tolerance_m3:raise ValueError('matched-footprint geometric/heat volume mismatch')
 record=dict(roi_id=ref.roi_id,bounds_m=list(ref.bounds_m),coordinate_frame=ref.frame,area_m2=float(ref.mask.sum()*ref.dx_m**2),dx_m=ref.dx_m,physical_time_s=physical_time_s,geometric_melt_m3=volume,heat_implied_melt_m3=heat_implied_volume_m3,geometry_heat_residual_m3=residual,absolute_ledger_tolerance_m3=absolute_ledger_tolerance_m3,mask_sha256=array_sha(ref.mask),initial_bed_sha256=array_sha(ref.bed_m),current_bed_sha256=array_sha(bed),provenance=ref.provenance,scope='fixed-XY ice-floor footprint; not vertical-cliff lifetime or field validation')
 return validate_record(record,mask=ref.mask,initial=ref.bed_m,current=bed)

def compare(candidate,reference,*,relative_reporting_min_volume_m3):
 validate_record(candidate);validate_record(reference)
 if not number(relative_reporting_min_volume_m3) or relative_reporting_min_volume_m3<=0:raise ValueError('predeclared positive relative reporting floor required')
 for k in ['roi_id','bounds_m','coordinate_frame','physical_time_s']:
  if candidate[k]!=reference[k]:raise ValueError('unmatched '+k)
 if abs(candidate['area_m2']-reference['area_m2'])>1e-12*max(1,reference['area_m2']):raise ValueError('unmatched physical footprint area')
 a=candidate['geometric_melt_m3'];b=reference['geometric_melt_m3']
 if not all(number(x) and x>=0 for x in [a,b]):raise ValueError('invalid melt volume')
 return dict(signed_absolute_error_m3=a-b,relative_error=(a/b-1 if b>=relative_reporting_min_volume_m3 else None),relative_reporting_min_volume_m3=relative_reporting_min_volume_m3,reference_is_observation=False,interpretation='matched-support model comparison only; source independence and physical validity are not certified')

def save_bundle(path,sim,ref,record,*,run_metadata):
 """Write float64 geometry+mask and manifest. Existing outputs never overwritten."""
 path=Path(path);data=path.with_suffix('.npz')
 if path.exists() or data.exists():raise ValueError('versioned output already exists')
 required={'solver_version','solver_source_sha256','configuration','evidence_class'}
 if not isinstance(run_metadata,dict) or not required<=set(run_metadata):raise ValueError('run metadata incomplete')
 if run_metadata['evidence_class']!='conditional_numerical_run':raise ValueError('this exporter cannot label field evidence')
 sha=run_metadata['solver_source_sha256']
 if not valid_sha(sha):raise ValueError('source SHA256 required')
 if not isinstance(run_metadata['solver_version'],str) or not run_metadata['solver_version'].strip() or not isinstance(run_metadata['configuration'],dict) or not run_metadata['configuration']:raise ValueError('version and complete configuration required')
 bed=np.asarray(sim.zb,dtype=np.float64);validate_record(record,mask=ref.mask,initial=ref.bed_m,current=bed)
 config_sha=digest_bytes(canonical(run_metadata['configuration']));canonical(run_metadata)
 buffer=io.BytesIO();np.savez_compressed(buffer,initial_bed_m=ref.bed_m,current_bed_m=bed,mask=ref.mask);raw=buffer.getvalue()
 manifest=dict(schema_version=1,data_file=data.name,data_sha256=digest_bytes(raw),record=record,run_metadata=run_metadata,configuration_sha256=config_sha,authentication='byte integrity only; not external run authentication')
 manifest_bytes=canonical(manifest)
 path.parent.mkdir(parents=True,exist_ok=True)
 def publish_no_overwrite(target,content):
  with tempfile.NamedTemporaryFile(mode='wb',dir=path.parent,delete=False) as f:
   tmp=Path(f.name);f.write(content);f.flush();os.fsync(f.fileno())
  try:os.link(tmp,target)
  finally:tmp.unlink()
 publish_no_overwrite(data,raw);publish_no_overwrite(path,manifest_bytes)
 return manifest

def load_bundle(path):
 path=Path(path);m=strict_json(path.read_text())
 if not isinstance(m,dict):raise ValueError('manifest object required')
 name=m.get('data_file')
 if m.get('schema_version')!=1 or not isinstance(name,str) or Path(name).name!=name:raise ValueError('untrusted/legacy bundle schema')
 data=path.parent/name
 if digest_bytes(data.read_bytes())!=m.get('data_sha256'):raise ValueError('bundle byte hash mismatch')
 meta=m.get('run_metadata')
 if not isinstance(meta,dict) or meta.get('evidence_class')!='conditional_numerical_run' or not valid_sha(meta.get('solver_source_sha256')) or not isinstance(meta.get('solver_version'),str) or not meta['solver_version'].strip():raise ValueError('invalid model provenance/evidence class')
 canonical(m)
 validate_record(m.get('record'))
 with zipfile.ZipFile(data) as archive:
  if sum(x.file_size for x in archive.infolist())>256*1024**2:raise ValueError('bundle exceeds supported uncompressed size')
 with np.load(data,allow_pickle=False) as z:
  if set(z.files)!={'initial_bed_m','current_bed_m','mask'}:raise ValueError('unexpected bundle arrays')
  a=z['initial_bed_m'];b=z['current_bed_m'];mask=z['mask'];r=m['record']
  if a.dtype!=np.float64 or b.dtype!=np.float64 or mask.dtype!=np.bool_ or a.shape!=b.shape or a.shape!=mask.shape:raise ValueError('float64 matching arrays required')
  if not np.isfinite(a).all() or not np.isfinite(b).all() or mask.ndim!=2 or mask.shape[0]<3 or mask.shape[1]<3:raise ValueError('invalid array contents')
  if mask[0].any() or mask[-1].any() or mask[:,0].any() or mask[:,-1].any():raise ValueError('ghost cells in ROI')
  if not number(r.get('dx_m')) or r['dx_m']<=0 or not number(r.get('physical_time_s')) or r['physical_time_s']<0:raise ValueError('invalid geometry/time metadata')
  if not mask.any() or not number(r.get('area_m2')) or abs(mask.sum()*r['dx_m']**2-r['area_m2'])>1e-12:raise ValueError('mask/area mismatch')
  if (b[mask]-a[mask]>1e-12).any():raise ValueError('unsupported ice growth')
  if array_sha(a)!=r['initial_bed_sha256'] or array_sha(b)!=r['current_bed_sha256'] or array_sha(mask)!=r['mask_sha256']:raise ValueError('record/array mismatch')
  volume=float(((a-b)[mask]).sum()*r['dx_m']**2)
  if abs(volume-r['geometric_melt_m3'])>1e-12*max(1,abs(volume)):raise ValueError('volume arithmetic mismatch')
  validate_record(r,mask=mask,initial=a,current=b)
 if digest_bytes(canonical(m['run_metadata']['configuration']))!=m['configuration_sha256']:raise ValueError('configuration hash mismatch')
 return m
