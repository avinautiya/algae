"""Recompute matched central-control melt comparisons from existing M3 exports.

Run from the packet folder:
 python3 -m melt.repair_eval EVAL_DIR production_spec.json --output matched.json
No CFD is launched. Missing raw records produce unavailable rows, never scores.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
from .observation import profile_melt,reduced_axial_roi,latent_energy

PIN='f3602e35e1595bcdbee4208a354a34944f1ae974'


def validate_output_path(eval_dir,spec_path,output_path):
    """Keep every output outside the raw-input tree and away from the spec."""
    output=Path(output_path).resolve();root=Path(eval_dir).resolve()
    if output==Path(spec_path).resolve():raise ValueError('output cannot overwrite spec')
    if output==root or root in output.parents:
        raise ValueError('output must be outside eval_dir; raw input files are immutable')
    if output.exists():
        spec=Path(spec_path)
        if spec.exists() and output.samefile(spec):raise ValueError('output aliases spec')
        if root.exists():
            for item in root.rglob('*'):
                if item.is_file() and output.samefile(item):raise ValueError('output aliases a raw input file')


def identity_check(raw,case):
    expected=dict(dx=case['dx'],layout=case['layout'],removal_day=case['removal_day'],
                  days=case['days'],dtm=case['dtm_h']*3600.)
    for field,value in expected.items():
        if field not in raw or raw[field]!=value:
            raise ValueError(case['name']+': raw '+field+' differs from frozen case or is missing')
        if isinstance(value,(int,float)) and (isinstance(raw[field],bool) or not math.isfinite(raw[field])):
            raise ValueError(case['name']+': nonnumeric/nonfinite raw '+field)


def repair(eval_dir,spec):
    root=Path(eval_dir);rows=[]
    for case in spec['evaluation_runs']:
        name=case['name'];path=root/(name+'_medium.json')
        row=dict(run=name,layout=case['layout'],status='unavailable',comparisons=[])
        if not path.exists():row['reason']='raw medium CFD record not supplied';rows.append(row);continue
        raw=json.loads(path.read_text())
        if raw.get('status')!='complete':row['reason']='CFD run not marked complete';rows.append(row);continue
        identity_check(raw,case)
        # The pinned campaign uses the same physical reach 3..5m and toe width .5m.
        # Do not silently reuse these coordinates for arbitrary future geometry.
        if case['dx'] not in (.025,.05,.1,.0125):raise ValueError('unsupported grid; supply a reviewed geometry adapter')
        if raw.get('records') and len(raw['records'][0].get('toe_floor_x',[]))!=round(2./case['dx']):
            raise ValueError(name+': exported reach profile does not match pinned2m geometry')
        try:
            control=profile_melt(raw['records'],dx_m=case['dx'],toe_width_m=.5,
                                 reach_start_m=3.,roi_start_m=3.5,roi_end_m=4.5)
            reach=profile_melt(raw['records'],dx_m=case['dx'],toe_width_m=.5,
                               reach_start_m=3.,roi_start_m=3.,roi_end_m=5.)
        except (ValueError,KeyError) as exc:
            row['reason']='missing/invalid CFD profiles: '+str(exc);rows.append(row);continue
        if control['first_day']!=0. or abs(control['last_day']-case['days'])>1e-9:
            row['reason']='CFD profile horizon differs from frozen prediction horizon';rows.append(row);continue
        initial_floor=spec.get('geometry',{}).get('z_toe_floor',-.08)
        if any(abs(z-initial_floor)>5e-6+1e-12 for z in raw['records'][0]['toe_floor_x']):
            row['reason']='initial toe profile differs from pinned initial floor';rows.append(row);continue
        row.update(status='reconstructed',cfd_control=control,cfd_full_reach=reach,
                   cfd_open_reach_melt_m3=reach['volume_m3']-control['volume_m3'],
                   raw_input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   provenance_status='metadata_identity_checked_not_run_authenticated')
        reported=raw.get('cum',{}).get('toe_reach')
        if reported is not None:
            difference=reach['volume_m3']-reported
            row['reported_cumulative_consistency']=dict(profile_minus_reported_m3=difference,
                exceeds_profile_rounding_bound=abs(difference)>reach['rounding_volume_bound_m3'],
                warning='profile lowering and nominal heat-based cumulative melt differ if ice-floor clipping occurs; investigate any excess')
        for model in spec.get('reduced_predictions',{}).get(name,[]):
            try:
                matched=reduced_axial_roi(model['toe_melt_m3'],method=model['method'],layout=case['layout'])
            except ValueError as exc:
                row['comparisons'].append(dict(method=model['method'],status='unsupported',reason=str(exc)));continue
            observed=control['volume_m3'];pred=matched['volume_m3'];error=pred-observed
            row['comparisons'].append(dict(method=model['method'],status='scope_matched',
                reduced_control=matched,cfd_control_m3=observed,error_m3=error,
                absolute_error_m3=abs(error),relative_error=error/observed if observed>control['rounding_volume_bound_m3'] else None,
                cfd_profile_rounding_bound_m3=control['rounding_volume_bound_m3'],
                reduced_latent_energy_J=latent_energy(pred),
                warning='common observation support only; shared closure and idealized reduced geometry remain'))
        rows.append(row)
    return dict(snapshot_commit=PIN,observation='central1m x full0.5m toe width, x3.5..4.5m',
                results=rows,scope='corrected footprint comparison; no new CFD or external validation',
                warning='Do not call unavailable rows zero error or aggregate shared trajectories as independent',
                provenance_limit='Filename and matching dx/layout/removal/day/dtm are not authentication. Raw runs omit Q,T and complete solver/config identity; production hash manifest and original execution provenance remain required.')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('eval_dir',type=Path);p.add_argument('spec',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    validate_output_path(a.eval_dir,a.spec,a.output)
    result=repair(a.eval_dir,json.loads(a.spec.read_text()))
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'output':str(a.output),'reconstructed':sum(r['status']=='reconstructed' for r in result['results']),
                      'unavailable':sum(r['status']=='unavailable' for r in result['results'])}))


if __name__=='__main__':main()
