"""Fixed-model future batch, with no new training and no camera outcome access.

Roll each already frozen physical column continuously from its original Dec2
initialization through the last source-eligible future endpoint. Never restart
at a test boundary, refit diffusivity/capacity/gains, or choose a new primary.
"""
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
from cold_ice_column import simulate
from piramide_experiment import process
from temperature_model import interpolate

ROOT=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def source(value):
    p=Path(value);return p if p.is_absolute() else ROOT/p
def load(value):return json.loads(source(value).read_text())


def interval_indices(labels,start,end):
    a,b=dt.datetime.fromisoformat(start),dt.datetime.fromisoformat(end)
    if not a<b or a not in labels or b not in labels:raise ValueError('exact increasing native thermal endpoints required')
    return labels.index(a),labels.index(b)


def sum_interval(rows,labels,start,end,C,L):
    a,b=interval_indices(labels,start,end)
    mass=math.fsum(r['melt_mmwe'] for r in rows[a:b])
    return dict(start=start,end=end,raw_melt_mmwe=mass,raw_latent_heat_J_m2=mass*L,
                normalized_basal_heat_K_m=mass*L/C,source_native_intervals=b-a)


def prepare(spec_path,output):
    if output.exists():raise ValueError('new batch-freeze directory required')
    c=json.loads(spec_path.read_text());v3f=load(c['source_v3_freeze']);v3=load(c['source_v3_normalized_predictions'])
    binding=load(c['endpoint_binding'])
    if binding['feature_coverage_json_sha256']!=sha(source(c['feature_coverage_audit'])):
        raise ValueError('endpoint binding/feature coverage audit mismatch')
    if binding['future_target_values_exposed'] is not False:raise ValueError('future target values must remain unopened')
    if sha(source(c['source_v3_freeze']))!=v3['freeze_sha256']:raise ValueError('source version3 prediction/freeze mismatch')
    eligible=[dict(start=r['start'],end=r['end']) for r in binding['intervals'] if r['source_eligible']]
    excluded=[dict(start=r['start'],end=r['end'],reason='source coverage failed before outcome access') for r in binding['intervals'] if not r['source_eligible']]
    if c['future_intervals']!=eligible or c['excluded_intervals']!=excluded:raise ValueError('batch must retain all and only source-eligible intervals')
    if len(eligible)!=3 or eligible[0]['start']!=v3['end_local'] or any(x['end']!=y['start'] for x,y in zip(eligible,eligible[1:])):
        raise ValueError('three consecutive future periods following source event required')
    if c['primary_metric']!='mean_absolute_period_ablation_error_mmwe':raise ValueError('fixed multi-period MAE required')
    if c['new_training_allowed'] is not False or c['new_model_selection_allowed'] is not False:raise ValueError('no retraining or model selection permitted')
    if c['primary_candidate']!='forward__middle_gap_20mm__left__proxy':raise ValueError('original primary remains fixed')
    if binding['selected_labels_unique'] is not True or binding['selected_labels_strictly_increasing'] is not True:
        raise ValueError('source endpoint collisions or ordering failure')
    # Copy physical family from the immutable already tested source. Config
    # contains lineage, not a replacement set of adjustable physical parameters.
    v1f=load(v3f['config']['v1_freeze'])
    times,_=process(source(v1f['config']['process_csv']),source(v1f['config']['process_manifest']))
    for row in eligible:interval_indices(times,row['start'],row['end'])
    if c['camera_evaluation']['intervals']!=[dict(tower=1,role='evaluation',**r) for r in eligible]:
        raise ValueError('camera receipt must bind only these future evaluation endpoints')
    source_names=['batch_phase_experiment.py','cold_ice_column.py','temperature_model.py','piramide_experiment.py',
                  'thermal_memory/solver.py','FUTURE_BATCH_PROTOCOL.md','FUTURE_BATCH_SPEC.json']
    inputs={k:sha(source(c[k])) for k in ('source_v3_freeze','source_v3_normalized_predictions',
              'source_v3_fixed_gain_predictions','endpoint_binding','feature_coverage_audit')}
    f=dict(config=c,source_hashes={n:sha(ROOT/n) for n in source_names},input_hashes=inputs,
           source_v3_process_hashes=v3f['v1_process_input_hashes'],
           physical_config=v3f['config'],camera_evaluation=c['camera_evaluation'],
           future_target_values_opened=False,future_model_predictions_computed=False,
           previous_single_event_result_known=True,
           primary_hypothesis='unchanged primary versus unchanged baselines on all three newly frozen future periods')
    output.mkdir(parents=True);(output/'FREEZE.json').write_text(json.dumps(f,indent=2)+'\n')
    return f


def run(directory,published_commit):
    target=directory/'NORMALIZED_HEAT_PREDICTIONS.json'
    if target.exists():raise ValueError('refusing to overwrite future-batch predictions')
    if len(published_commit)!=40 or any(x not in '0123456789abcdef' for x in published_commit):raise ValueError('published pre-execution commit required')
    f=json.loads((directory/'FREEZE.json').read_text());c=f['config'];physical=f['physical_config']
    for name,expected in f['source_hashes'].items():
        if sha(ROOT/name)!=expected:raise ValueError('frozen batch model/protocol changed')
    for name,expected in f['input_hashes'].items():
        if sha(source(c[name]))!=expected:raise ValueError('frozen batch source lineage changed')
    v3=load(c['source_v3_normalized_predictions']);v1=load(physical['v1_freeze'])
    for name,expected in f['source_v3_process_hashes'].items():
        if sha(source(v1['config'][name]))!=expected:raise ValueError('original source/process input changed')
    times,values=process(source(v1['config']['process_csv']),source(v1['config']['process_manifest']))
    initial=dt.datetime.fromisoformat(v3['start_local']);finish=dt.datetime.fromisoformat(c['future_intervals'][-1]['end'])
    ids=[i for i,t in enumerate(times) if initial<=t<=finish];labels=[times[i] for i in ids];data=[values[i] for i in ids]
    if labels[0]!=initial or labels[-1]!=finish:raise ValueError('continuous exact source support required')
    if len(v3['results'])!=36:raise ValueError('all36 physical family members required')
    results=[];L=physical['latent_heat_J_kg']
    for old in v3['results']:
        if old['status']!='predicted':raise ValueError('source family contains a failed member; no silent exclusion')
        g=old['geometry'];top=g['sensor_depths_m'][0];depth=g['ice_depth_m']-top
        zz=[z-top for z in g['sensor_depths_m']];n=physical['debris_cells']
        Td=interpolate(zz+[depth],data[0]+[0.],[(j+.5)*depth/n for j in range(n)])
        upper=[r[0] for r in (data[:-1] if old['hold']=='left' else data[1:])]
        solved=simulate(debris_conductivity_W_mK=old['debris_conductivity_W_mK'],
            debris_heat_capacity_J_m3K=old['physicalcapacity_J_m3K'],debris_depth_m=depth,
            initial_debris_temperature_C=Td,initial_ice_temperature_C=[physical['cold_ice_initial_temperature_C']]*physical['ice_cells'],
            surface_temperature_C=upper,interval_seconds=1800.,dt_seconds=physical['dt_seconds'],
            ice_depth_m=physical['ice_depth_m'],ice_conductivity_W_mK=physical['ice_conductivity_W_mK'],
            ice_heat_capacity_J_m3K=physical['ice_heat_capacity_J_m3K'],latent_heat_J_kg=L)
        intervals=solved['intervals'];C=old['physicalcapacity_J_m3K']
        past_train=sum_interval(intervals,labels,v3['start_local'],v3['split_local'],C,L)
        past_eval=sum_interval(intervals,labels,v3['split_local'],v3['end_local'],C,L)
        for key,actual in [('raw_training_melt_mmwe',past_train['raw_melt_mmwe']),('raw_evaluation_melt_mmwe',past_eval['raw_melt_mmwe'])]:
            if not math.isclose(actual,old[key],rel_tol=1e-11,abs_tol=1e-10):raise ValueError('continuous source prefix no longer reproduces version3: '+key)
        blocks=[sum_interval(intervals,labels,r['start'],r['end'],C,L) for r in c['future_intervals']]
        prefix='forward' if old['method']=='training_temperature_forward_fit' else 'derivative_OLS'
        ident=prefix+'__'+g['id']+'__'+old['hold']+'__'+old['capacity_id']
        results.append(dict(model_id=ident,method=old['method'],geometry=g,hold=old['hold'],capacity_id=old['capacity_id'],
            physicalcapacity_J_m3K=C,alpha_m2_s=old['alpha_m2_s'],debris_conductivity_W_mK=old['debris_conductivity_W_mK'],
            status='predicted',future_blocks=blocks,source_version3_prefix_reproduced=True,
            original_training_normalized_basal_heat_K_m=past_train['normalized_basal_heat_K_m'],
            prior_known_evaluation_normalized_basal_heat_K_m=past_eval['normalized_basal_heat_K_m'],
            cumulative_raw_melt_mmwe=solved['melt_mmwe'],cumulative_raw_ice_equivalent_retreat_m=solved['melt_mmwe']/917.,
            physical_energy_ledger=dict(initial_sensible_storage_J_m2=solved['initial_sensible_storage_J_m2'],
                final_sensible_storage_J_m2=solved['final_sensible_storage_J_m2'],surface_heat_J_m2=solved['surface_heat_J_m2'],
                deep_heat_out_J_m2=solved['deep_heat_out_J_m2'],latent_heat_J_m2=solved['latent_heat_J_m2'],
                energy_residual_J_m2=solved['energy_residual_J_m2'],relative_energy_residual=solved['relative_energy_residual'],
                max_step_energy_residual_J_m2=solved['max_step_energy_residual_J_m2']),
            complementarity_max_product=solved['complementarity_max_product_C_kg_m2'],
            minimum_interface_temperature_C=min(x['interface_temperature_C'] for x in solved['history']),
            cold_interface_substeps=sum(x['branch']=='cold_interface' for x in solved['history']),
            raw_component_ledger_does_not_apply_to_empirically_scaled_outputs=True,
            no_refitting_or_test_boundary_state_reset=True))
    result=dict(preexecution_commit=published_commit,freeze_sha256=sha(directory/'FREEZE.json'),
        continuous_model_initialization_local=initial.isoformat(),continuous_model_end_local=finish.isoformat(),
        future_intervals=c['future_intervals'],excluded_intervals=c['excluded_intervals'],results=results,
        primary_candidate=c['primary_candidate'],primary_metric=c['primary_metric'],
        source_fixed_gain_predictions_sha256=f['input_hashes']['source_v3_fixed_gain_predictions'],
        gains_refitted=False,primary_changed=False,future_camera_values_opened=False,
        previous_single_event_failure_preserved=True,
        scope='Conditional observed-temperature-forced rolling batch on one tower; fixed geometry, no wet-pore enthalpy or ice advection')
    with target.open('x') as stream:stream.write(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return dict(output=str(target),sha256=sha(target),physical_models=len(results),future_periods=len(c['future_intervals']),future_camera_values_opened=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['freeze','predict']);p.add_argument('path',type=Path)
    p.add_argument('--out',type=Path);p.add_argument('--commit');a=p.parse_args()
    try:
        if a.mode=='freeze':
            if not a.out:p.error('--out NEW required')
            prepare(a.path,a.out);print('future batch freeze written; no future prediction or outcome read')
        else:print(json.dumps(run(a.path,a.commit or '')))
    except (ValueError,KeyError,TypeError) as error:p.error(str(error))
