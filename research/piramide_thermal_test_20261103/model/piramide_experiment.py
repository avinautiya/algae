"""Process-only, frozen conditional ice-column prediction campaign.

This module never opens camera MAT files or camera outcomes. It outputs native
normalized basal-heat traces for a separately sealed, training-only melt gain.
"""
import argparse,csv,datetime as dt,hashlib,json,math
from pathlib import Path
from inference import fit_grid,derivative_crh
from temperature_model import predict

ROOT=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stamp(text):
    t=dt.datetime.fromisoformat(text)
    if t.tzinfo is not None:raise ValueError('native clock is naive; do not assign UTC')
    return t


def process(path,manifest):
    record=json.loads(Path(manifest).read_text())
    entry=next(x for x in record['checks'] if x['family']=='debtemp' and x['tower']=='Tower1')
    if sha(path)!=entry['sha256']:raise ValueError('verified source process hash mismatch')
    with Path(path).open() as f:
        reader=csv.DictReader(f)
        if reader.fieldnames!=['SOURCE_LOCAL_LABEL','DT1','DT2','DT3','DT4','DT5']:
            raise ValueError('native process-only temperature schema required')
        rows=list(reader)
    times=[stamp(x['SOURCE_LOCAL_LABEL']) for x in rows]
    values=[[float(x['DT'+str(j)]) for j in range(1,6)] for x in rows]
    if len(times)!=entry['rows'] or any((b-a).total_seconds()!=1800 for a,b in zip(times,times[1:])):
        raise ValueError('complete native 30-minute time support required')
    if not all(math.isfinite(v) for row in values for v in row):raise ValueError('missing/nonfinite temperature refused')
    return times,values


def check_config(c):
    if c.get('tower')!='Tower1':raise ValueError('primary Tower1 contract required')
    if c.get('clock')!='same_source_naive_local':raise ValueError('common source clock contract required')
    if c.get('candidate_claim')!='conditional_melt_predictor_with_train_only_gain':raise ValueError('conditional claim contract required')
    if not c.get('physical_assumptions'):raise ValueError('explicit physical assumptions required')
    if c.get('camera_outcomes_opened_before_freeze') is not False:raise ValueError('pre-outcome freeze required')
    a,b,z=[stamp(c[k]) for k in ('start','split','end')]
    if not a<b<z:raise ValueError('chronological disjoint periods required')
    if type(c['cells']) is not int or c['cells']<4:raise ValueError('integer cells>=4 required')
    if type(c['dt_seconds']) not in (int,float) or not math.isfinite(c['dt_seconds']) or c['dt_seconds']<=0:
        raise ValueError('finite positive dt required')
    grid=c['alpha_grid_m2_s']
    if not grid or any(type(v) not in (int,float) or not math.isfinite(v) or v<=0 for v in grid) or any(v<=u for u,v in zip(grid,grid[1:])):
        raise ValueError('explicit increasing positive alpha grid required')
    ids=[];nominal=[0,.125,.25,.375,.5]
    for scene in c['geometries']:
        ids.append(scene['id']);d=scene['ice_depth_m'];zz=scene['sensor_depths_m']
        if len(zz)!=5 or not all(type(v) in (int,float) and math.isfinite(v) for v in zz+[d]):raise ValueError('finite geometry required')
        if any(abs(x-y)>.020000000001 for x,y in zip(zz,nominal)) or abs(d-.5)>.020000000001:
            raise ValueError('geometry exceeds published ±2cm bounds')
        if any(v<=u for u,v in zip(zz,zz[1:])) or zz[-1]>=d:raise ValueError('DT5 must be interior above ice with positive gap')
        if not isinstance(scene.get('source'),str) or not scene['source'].strip():raise ValueError('geometry source/assumption required')
    if not ids or len(set(ids))!=len(ids):raise ValueError('nonempty uniquely named geometry family required')
    if c['holds']!=['left','right']:raise ValueError('both native reconstruction scenarios required; do not choose winner')
    if c.get('primary_geometry') not in ids or c.get('primary_hold') not in c['holds']:
        raise ValueError('primary geometry/hold must be fixed before targets')
    camera=c['camera_evaluation']
    expected_intervals=[dict(tower=1,role='train',start=c['start'],end=c['split']),dict(tower=1,role='evaluation',start=c['split'],end=c['end'])]
    if camera.get('intervals')!=expected_intervals:raise ValueError('camera train/eval endpoints disagree')
    if camera.get('approved_operator')!='exposed_stake_distance_endpoint_difference_m_ice_times_917_to_mmwe':
        raise ValueError('source raw distance operator required')
    if camera.get('clock_interpretation')!='same_source_naive_local_no_UTC_assumption':raise ValueError('native camera clock contract required')
    if camera.get('camera_manifest_sha256')!=sha(c['camera_metadata']):raise ValueError('sealed camera manifest hash mismatch')
    return a,b,z


def prepare(config_path,output):
    if output.exists():raise ValueError('new freeze directory required')
    raw=Path(config_path).read_bytes();c=json.loads(raw);a,b,z=check_config(c)
    times,values=process(c['process_csv'],c['process_manifest'])
    if not all(t in times for t in (a,b,z)):raise ValueError('exact split endpoints require native temperature support')
    for name in ('camera_times_only','camera_metadata'):
        if not Path(c[name]).is_file():raise ValueError('sealed camera clock/metadata file required')
    # Calendar endpoints must be verified by the separately audited timestamp-only
    # source. No native camera displacement/ablation field is loaded here.
    camera_times=json.loads(Path(c['camera_times_only']).read_text())
    labels=camera_times if isinstance(camera_times,list) else camera_times.get('source_local_times',camera_times.get('timestamps',[]))
    labels={stamp(t) for t in labels}
    if not all(t in labels for t in (a,b,z)):raise ValueError('exact camera timestamp support missing; no endpoint interpolation')
    names=['piramide_experiment.py','temperature_model.py','inference.py','thermal_memory/solver.py',
           'PROTOCOL_PIRAMIDE.md','PIRAMIDE_EXPERIMENT_SPEC.json']
    source_hashes={n:sha(ROOT/n) for n in names}
    input_hashes={k:sha(c[k]) for k in ('process_csv','process_manifest','camera_times_only','camera_metadata')}
    freeze=dict(config=c,config_sha256=hashlib.sha256(raw).hexdigest(),source_hashes=source_hashes,
                input_hashes=input_hashes,process_camera_targets_read=False,training_fit_run=False,
                heldout_predictions_computed=False,clock='same-source naive local; no inferred UTC',
                endpoint_operator='exact source timestamp endpoint differences; camera units still independently gated',
                camera_evaluation=c['camera_evaluation'],
                independent_model_novelty_claim=False,professional_superiority_claim=False)
    output.mkdir(parents=True);(output/'FREEZE.json').write_text(json.dumps(freeze,indent=2)+'\n')
    return freeze


def predict_frozen(directory,preexecution_commit):
    out=directory/'NORMALIZED_HEAT_PREDICTIONS.json'
    if out.exists():raise ValueError('refusing to overwrite predictions')
    if len(preexecution_commit)!=40 or any(v not in '0123456789abcdef' for v in preexecution_commit):
        raise ValueError('published pre-execution commit SHA required')
    frozen=json.loads((directory/'FREEZE.json').read_text());c=frozen['config'];start,split,end=check_config(c)
    for n,h in frozen['source_hashes'].items():
        if sha(ROOT/n)!=h:raise ValueError('model source changed after freeze')
    for k,h in frozen['input_hashes'].items():
        if sha(c[k])!=h:raise ValueError('source/clock input changed after freeze')
    times,values=process(c['process_csv'],c['process_manifest'])
    indices=[i for i,t in enumerate(times) if start<=t<=end]
    selected_times=[times[i] for i in indices];v=[values[i] for i in indices]
    seconds=[(t-start).total_seconds() for t in selected_times]
    train_n=selected_times.index(split)
    train_times=seconds[:train_n];train_profiles=v[:train_n]
    fit_end=train_times[-1];evaluation_start=seconds[train_n]
    results=[]
    for geometry in c['geometries']:
        absolute=geometry['sensor_depths_m'];top=absolute[0]
        sensors=[z-top for z in absolute];depth=geometry['ice_depth_m']-top
        for hold in c['holds']:
            def solve(alpha,rows):
                return predict(alpha,domain_depth_m=depth,
                    initial_depths_m=sensors+[depth],initial_temperature_C=rows[0]+[0.],
                    upper_temperature_C=[r[0] for r in rows],lower_temperature_C=[0.]*len(rows),
                    interval_seconds=1800.,sample_depths_m=sensors[1:],cells=c['cells'],
                    dt_seconds=c['dt_seconds'],hold=hold)
            fitted=fit_grid(train_times,[row[1:] for row in train_profiles],sensor_ids=['DT2','DT3','DT4','DT5'],
                boundary_sensor_ids=['DT1'],train_end_s=fit_end,evaluation_start_s=evaluation_start,
                alpha_grid_m2_s=c['alpha_grid_m2_s'],predict_training=lambda alpha:solve(alpha,train_profiles)['predicted_temperature_C'],
                sensor_absolute_bound_C=c['sensor_absolute_bound_C'])
            methods=[('training_temperature_forward_fit',fitted['selected_alpha_m2_s'],fitted)]
            try:
                derivative=derivative_crh(train_times,train_profiles,sensors,train_end_s=fit_end,evaluation_start_s=evaluation_start)
                methods.append(('training_derivative_OLS_NOT_source_NTLS',derivative['alpha_m2_s'],derivative))
            except ValueError as error:
                results.append(dict(geometry=geometry,hold=hold,method='training_derivative_OLS_NOT_source_NTLS',status='refused',reason=str(error)))
            for method,alpha,inference in methods:
                solved=solve(alpha,v);heat=solved['interval_basal_normalized_heat_K_m']
                trace=[dict(start_s=seconds[i],end_s=seconds[i+1],start_local=selected_times[i].isoformat(),
                            end_local=selected_times[i+1].isoformat(),normalized_basal_heat_K_m=h)
                       for i,h in enumerate(heat)]
                train_heat=math.fsum(heat[:train_n]);eval_heat=math.fsum(heat[train_n:])
                results.append(dict(geometry=geometry,hold=hold,method=method,status='predicted',alpha_m2_s=alpha,
                    inference=inference,normalized_heat_trace=trace,
                    training_normalized_basal_heat_K_m=train_heat,evaluation_normalized_basal_heat_K_m=eval_heat,
                    evaluation_training_heat_ratio=(eval_heat/train_heat if train_heat>0 else None),
                    negative_basal_heat_intervals=sum(h<0 for h in heat),
                    negative_basal_flux_substeps=solved['negative_basal_flux_substeps'],
                    no_refreezing_melt_conversion_eligible=all(h>=0 for h in heat) and train_heat>0 and solved['negative_basal_flux_substeps']==0,
                    normalized_energy_residual_K_m=solved['normalized_energy_residual_K_m'],
                    predicted_internal_temperature_C=solved['predicted_temperature_C'],
                    measured_evaluation_internal_temperature_scored=False))
    result=dict(preexecution_commit=preexecution_commit,freeze_sha256=sha(directory/'FREEZE.json'),
                start_local=start.isoformat(),split_local=split.isoformat(),end_local=end.isoformat(),
                relative_train_end_s=(split-start).total_seconds(),relative_evaluation_start_s=(split-start).total_seconds(),
                process_native_cadence_seconds=1800,results=results,
                primary_geometry=c['primary_geometry'],primary_hold=c['primary_hold'],
                camera_outcome_fields_loaded=False,calibration_melt_gain_fitted=False,
                scope='Conditional ice-to-buried-thermistor geometry family. Normalized heat is not yet absolute melt.',
                assumptions=c['physical_assumptions'],model_novelty_claim=False,field_accuracy_established=False)
    with out.open('x') as f:f.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return dict(output=str(out),sha256=sha(out),scenario_methods=len(results),camera_outcomes_loaded=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['freeze','predict']);p.add_argument('path',type=Path)
    p.add_argument('--out',type=Path);p.add_argument('--commit');a=p.parse_args()
    try:
        if a.mode=='freeze':
            if not a.out:p.error('--out NEW_DIR required')
            prepare(a.path,a.out);print('freeze written; no model fit or camera outcomes loaded')
        else:print(json.dumps(predict_frozen(a.path,a.commit or '')))
    except (ValueError,KeyError,TypeError,StopIteration) as error:p.error(str(error))


if __name__=='__main__':main()
