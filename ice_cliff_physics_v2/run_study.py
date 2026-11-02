"""Versioned synthetic study with source/input freezes and atomic case records."""
import argparse,datetime,hashlib,json,math,os,sys,tempfile,time,traceback
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from model import forecast


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def source_files():
    paths=[ROOT/'model.py',ROOT/'face_thermal.py',Path(__file__),ROOT/'run_boundary_study.py',ROOT/'mechanism/wet_resistance.py']
    for folder in ['thermal','transfer','face_geometry','energy_ports','vendor']:
        paths.extend(p for p in (ROOT/folder).glob('*.py') if not p.name.startswith('test_'))
    return {p.relative_to(ROOT).as_posix():digest(p) for p in sorted(set(paths))}

def atomic(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,delete=False) as f:
        name=f.name;json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(name,path)

def freeze():
    path=ROOT/'study/FREEZE.json'
    if path.exists():raise ValueError('freeze exists; retain it and create a versioned new study for any amendment')
    record=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scope='synthetic verification and approximation comparison; not blind field validation',
        coupled_outcomes_run=False,component_tests_and_analytic_limits_already_seen=True,
        cases_sha256=digest(ROOT/'study/COUPLED_CASES.json'),spec_sha256=digest(ROOT/'study/SPEC.md'),
        wet_boundary_cases_sha256=digest(ROOT/'mechanism/wet_boundary_cases.json'),sources=source_files())
    atomic(path,record);print(json.dumps(record,indent=2),flush=True)

def summarize(result,config):
    series=result['series'];physics=config['physics'];geometry=config['geometry']
    covered_area=(config['upgrades']['fixed_cover_interval_m'][1]-config['upgrades']['fixed_cover_interval_m'][0])*geometry['extrusion_width_m'] if 'fixed_cover_interval_m' in config.get('upgrades',{}) else None
    fixed=(covered_area is None or all(abs(row['port_states']['face_debris']['area_m2']-covered_area)<1e-12 for row in series))
    sources=max(1.,abs(result['totals']['atmospheric_J'])+result['totals']['water_in_J']+abs(result['totals']['remap_enthalpy_J']))
    normalized=max(abs(r['global_energy_residual_J'])/r['cumulative_energy_accounting_magnitude_J'] for r in series)
    if normalized>1e-8 or abs(result['checks']['geometry_phase_volume_residual_m3'])>1e-9 or not fixed:
        raise RuntimeError('frozen numerical or matched-footprint acceptance criterion failed')
    available=[i for i,r in enumerate(series) if r['face_profile'] is not None]
    selected=sorted(set(min(available,key=lambda i:abs(i-target)) for target in [0,len(series)//4,len(series)//2,3*len(series)//4,len(series)-1])) if available else [0,len(series)-1]
    snapshots=[dict(time_s=series[i]['time_s'],routing_geometry=series[i]['geometry'],
                    face_profile=series[i]['face_profile'],undercut=series[i]['undercut_diagnostic'],
                    debris_columns=series[i].get('debris_columns_post_geometry')) for i in selected]
    trace=[]
    for r in series:
        trace.append(dict(time_s=r['time_s'],covered_melt_m3=r['covered_face_melt_m3'],
            face_melt_m3=r['face_melt_m3'],floor_melt_m3=r['floor_melt_m3'],
            toe_flow_m3_s=r['toe_flow_m3_s'],toe_fraction=r['toe_fraction'],
            energy_residual_J=r['global_energy_residual_J'],
            profile_bands=r['profile_height_band_count'],
            undercut_m=0. if r['undercut_diagnostic'] is None else r['undercut_diagnostic']['maximum_horizontal_undercut_m']))
    return dict(totals=result['totals'],checks=result['checks'],fixed_covered_area_m2=covered_area,
                matched_covered_footprint_pass=fixed,max_time_series_relative_energy_residual=normalized,
                maximum_time_series_energy_residual_J=max(abs(r['global_energy_residual_J']) for r in series),
                time_steps=len(series),maximum_profile_bands=max((r['profile_height_band_count'] or 0) for r in series),
                snapshots=snapshots,trace=trace,scope=result['scope'])

def run(out,commit):
    if len(commit)!=40 or any(c not in '0123456789abcdef' for c in commit):raise ValueError('published pre-execution commit SHA required')
    record=json.loads((ROOT/'study/FREEZE.json').read_text())
    if record['sources']!=source_files() or record['cases_sha256']!=digest(ROOT/'study/COUPLED_CASES.json') or record['spec_sha256']!=digest(ROOT/'study/SPEC.md'):
        raise ValueError('source/spec/input changed after freeze; no silent reuse')
    cases=json.loads((ROOT/'study/COUPLED_CASES.json').read_text());summary=[]
    for case in cases:
        path=out/(case['id']+'.json');input_hash=hashlib.sha256(json.dumps(case,sort_keys=True,allow_nan=False).encode()).hexdigest()
        if path.exists():
            old=json.loads(path.read_text())
            if old['case_sha256']!=input_hash or old['freeze_sha256']!=digest(ROOT/'study/FREEZE.json'):
                raise ValueError('existing case belongs to a different frozen run')
            summary.append({k:old[k] for k in ['id','status','runtime_s']});continue
        start=time.perf_counter()
        result=dict(id=case['id'],case_sha256=input_hash,freeze_sha256=digest(ROOT/'study/FREEZE.json'),
                    freeze_commit=commit,configuration=case['configuration'],backend=case['backend'],
                    scope='synthetic coupled model comparison, not environmental observations')
        try:
            prediction=forecast(case['configuration']);result.update(status='completed',result=summarize(prediction,case['configuration']))
        except ValueError as exc:result.update(status='refused',reason=str(exc))
        except Exception as exc:result.update(status='failed',reason=str(exc),traceback=traceback.format_exc())
        result['runtime_s']=time.perf_counter()-start;atomic(path,result)
        summary.append({k:result[k] for k in ['id','status','runtime_s']})
        print(case['id'],result['status'],round(result['runtime_s'],3),flush=True)
        if result['status']!='completed':print(result['reason'],flush=True)
    atomic(out/'RUN_INDEX.json',dict(freeze_commit=commit,cases=summary,
        completed=sum(x['status']=='completed' for x in summary),refused=sum(x['status']=='refused' for x in summary),failed=sum(x['status']=='failed' for x in summary)))
    if any(x['status']=='failed' for x in summary):raise SystemExit(1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true');p.add_argument('--out',type=Path);p.add_argument('--freeze-commit');a=p.parse_args()
    if a.freeze:freeze()
    elif a.out and a.freeze_commit:run(a.out,a.freeze_commit)
    else:p.error('use --freeze, or --out DIRECTORY --freeze-commit SHA')
