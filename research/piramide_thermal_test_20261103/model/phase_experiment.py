"""Version2 cold-ice extension; preserves the published version1 phase failure.

Produces raw physical melt and normalized latent-energy traces without reading
camera outcomes. Any later training gain is an empirical output correction,
not a re-identification of physical heat capacity.
"""
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
from cold_ice_column import simulate
from piramide_experiment import process
from temperature_model import interpolate

ROOT=Path(__file__).resolve().parent
CAPACITY_KEYS={'dry':'dry_assumed_pore_air','proxy':'training_valid_probe_proxy_not_measured_depth_mean',
               'saturated':'fully_water_saturated_conditional'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def path(value):
    x=Path(value);return x if x.is_absolute() else ROOT/x


def prepare(spec_path,output):
    if output.exists():raise ValueError('new version2 freeze directory required')
    c=json.loads(spec_path.read_text());v1=json.loads(path(c['v1_freeze']).read_text())
    normalized=json.loads(path(c['v1_normalized_predictions']).read_text())
    if normalized['freeze_sha256']!=sha(path(c['v1_freeze'])):raise ValueError('version1 artifact/freeze mismatch')
    if c['alpha_source']!='unchanged_v1_training_only_temperature_fits':raise ValueError('explicit alpha reuse required')
    if c['cold_ice_initial_temperature_C']!=0.:raise ValueError('version2 primary initial ice is0C; changes require new spec')
    if c['primary_capacity_id']!='proxy':raise ValueError('source training-proxy primary required')
    if c['primary_geometry']!=v1['config']['primary_geometry'] or c['primary_hold']!=v1['config']['primary_hold']:
        raise ValueError('no outcome-guided primary geometry/reconstruction change')
    for key in ('start','split','end'):
        if c[key]!=v1['config'][key]:raise ValueError('same declared chronological windows required')
    if c['camera_evaluation']!=v1['camera_evaluation']:raise ValueError('camera operator/roles must remain unchanged')
    audit=json.loads(path(c['capacity_source_audit']).read_text())
    source_window=[dt.datetime.fromisoformat(x) for x in audit['training_interval_half_open_source_local']]
    if source_window!=[dt.datetime.fromisoformat(c[k]) for k in ('start','split')]:raise ValueError('capacity source must use TRAIN only')
    expected=[dict(id=k,capacity_J_m3K=audit['finite_sensible_capacity_scenes'][v]['C_J_m3K']) for k,v in CAPACITY_KEYS.items()]
    if c['capacities']!=expected:raise ValueError('capacity scenes differ from independent pre-outcome TRAIN audit')
    if c['geometries']!=v1['config']['geometries'] or c['holds']!=v1['config']['holds']:
        raise ValueError('full original geometry/hold family required')
    names=['phase_experiment.py','cold_ice_column.py','temperature_model.py','piramide_experiment.py',
           'PHASE_EXTENSION_PROTOCOL.md','PHASE_EXTENSION_SPEC.json','thermal_memory/solver.py']
    inputs={key:sha(path(c[key])) for key in ('v1_freeze','v1_normalized_predictions','capacity_source_audit')}
    for key,expected_sha in v1['input_hashes'].items():
        source=path(v1['config'][key])
        if sha(source)!=expected_sha:raise ValueError('original process/clock input changed')
    f=dict(config=c,source_hashes={n:sha(ROOT/n) for n in names},input_hashes=inputs,
           v1_process_input_hashes=v1['input_hashes'],camera_evaluation=c['camera_evaluation'],
           physical_extension='fixed-geometry cold ice complementarity; not debris pore enthalpy',
           camera_outcomes_read=False,new_field_model_computed=False)
    output.mkdir(parents=True);(output/'FREEZE.json').write_text(json.dumps(f,indent=2)+'\n')
    return f


def run(directory,published_commit):
    if len(published_commit)!=40 or any(x not in '0123456789abcdef' for x in published_commit):raise ValueError('published pre-execution commit required')
    target=directory/'NORMALIZED_HEAT_PREDICTIONS.json'
    if target.exists():raise ValueError('refusing to overwrite version2 predictions')
    f=json.loads((directory/'FREEZE.json').read_text());c=f['config']
    for name,expected in f['source_hashes'].items():
        if sha(ROOT/name)!=expected:raise ValueError('version2 frozen source changed')
    for name,expected in f['input_hashes'].items():
        if sha(path(c[name]))!=expected:raise ValueError('version2 provenance input changed')
    v1=json.loads(path(c['v1_freeze']).read_text());old=json.loads(path(c['v1_normalized_predictions']).read_text())
    for name,expected in f['v1_process_input_hashes'].items():
        if sha(path(v1['config'][name]))!=expected:raise ValueError('process source changed')
    times,values=process(path(v1['config']['process_csv']),path(v1['config']['process_manifest']))
    start,split,end=[dt.datetime.fromisoformat(c[k]) for k in ('start','split','end')]
    selected=[i for i,t in enumerate(times) if start<=t<=end];labels=[times[i] for i in selected];data=[values[i] for i in selected]
    if labels[0]!=start or labels[-1]!=end or split not in labels:raise ValueError('exact process endpoints required')
    seconds=[(t-start).total_seconds() for t in labels];split_index=labels.index(split);results=[]
    for original in old['results']:
        for capacity in c['capacities']:
            if original['status']!='predicted':
                results.append(dict(method=original['method'],geometry=original['geometry'],hold=original['hold'],
                                    capacity_id=capacity['id'],status='refused',reason='original alpha method refused'));continue
            geom=original['geometry'];top=geom['sensor_depths_m'][0]
            depths=[z-top for z in geom['sensor_depths_m']];depth=geom['ice_depth_m']-top
            cells=c['debris_cells'];centers=[(i+.5)*depth/cells for i in range(cells)]
            initial=interpolate(depths+[depth],data[0]+[0.],centers)
            C=capacity['capacity_J_m3K'];alpha=original['alpha_m2_s']
            upper=[row[0] for row in (data[:-1] if original['hold']=='left' else data[1:])]
            solved=simulate(debris_conductivity_W_mK=alpha*C,debris_heat_capacity_J_m3K=C,
                debris_depth_m=depth,initial_debris_temperature_C=initial,
                initial_ice_temperature_C=[c['cold_ice_initial_temperature_C']]*c['ice_cells'],
                surface_temperature_C=upper,interval_seconds=1800.,dt_seconds=c['dt_seconds'],
                ice_depth_m=c['ice_depth_m'],ice_conductivity_W_mK=c['ice_conductivity_W_mK'],
                ice_heat_capacity_J_m3K=c['ice_heat_capacity_J_m3K'],latent_heat_J_kg=c['latent_heat_J_kg'])
            trace=[]
            for i,row in enumerate(solved['intervals']):
                latent=row['melt_kg_m2']*c['latent_heat_J_kg'];heat=latent/C
                trace.append(dict(start_s=seconds[i],end_s=seconds[i+1],start_local=labels[i].isoformat(),
                    end_local=labels[i+1].isoformat(),normalized_basal_heat_K_m=heat,
                    raw_melt_mmwe=row['melt_mmwe'],raw_latent_heat_J_m2=latent))
            train=math.fsum(r['normalized_basal_heat_K_m'] for r in trace[:split_index])
            evaluation=math.fsum(r['normalized_basal_heat_K_m'] for r in trace[split_index:])
            results.append(dict(method=original['method'],geometry=geom,hold=original['hold'],
                capacity_id=capacity['id'],capacity_J_m3K=C,physicalcapacity_J_m3K=C,status='predicted',alpha_m2_s=alpha,
                debris_conductivity_W_mK=alpha*C,alpha_refitted=False,normalized_heat_trace=trace,
                training_normalized_basal_heat_K_m=train,evaluation_normalized_basal_heat_K_m=evaluation,
                evaluation_training_heat_ratio=evaluation/train if train>0 else None,
                raw_training_melt_mmwe=math.fsum(r['raw_melt_mmwe'] for r in trace[:split_index]),
                raw_evaluation_melt_mmwe=math.fsum(r['raw_melt_mmwe'] for r in trace[split_index:]),
                negative_basal_heat_intervals=0,negative_latent_melt_substeps=0,
                negative_basal_flux_substeps=sum(r['debris_interface_heat_W_m2']<0 for r in solved['history']),
                negative_basal_flux_is_sensible_cooling_supported_by_ice_domain=True,
                no_refreezing_melt_conversion_eligible=train>0,
                cold_interface_substeps=sum(r['branch']=='cold_interface' for r in solved['history']),
                minimum_interface_temperature_C=min(r['interface_temperature_C'] for r in solved['history']),
                interface_complementarity_max_product=solved['complementarity_max_product_C_kg_m2'],
                physical_energy_ledger=dict(initial_sensible_storage_J_m2=solved['initial_sensible_storage_J_m2'],
                    final_sensible_storage_J_m2=solved['final_sensible_storage_J_m2'],surface_heat_J_m2=solved['surface_heat_J_m2'],
                    deep_heat_out_J_m2=solved['deep_heat_out_J_m2'],latent_heat_J_m2=solved['latent_heat_J_m2'],
                    energy_residual_J_m2=solved['energy_residual_J_m2'],relative_energy_residual=solved['relative_energy_residual'],
                    max_step_energy_residual_J_m2=solved['max_step_energy_residual_J_m2']),
                ledger_applies_to_unscaled_raw_component_only=True,
                empirical_output_gain_is_not_physical_heat_capacity=True,
                wet_debris_phase_change_included=False,fixed_geometry_no_ice_thermal_advection=True))
    report=dict(preexecution_commit=published_commit,freeze_sha256=sha(directory/'FREEZE.json'),
        start_local=start.isoformat(),split_local=split.isoformat(),end_local=end.isoformat(),
        relative_train_end_s=seconds[split_index],relative_evaluation_start_s=seconds[split_index],
        process_native_cadence_seconds=1800,results=results,primary_geometry=c['primary_geometry'],
        primary_hold=c['primary_hold'],primary_capacity_id=c['primary_capacity_id'],
        camera_outcome_fields_loaded=False,calibration_melt_gain_fitted=False,
        raw_component_uses_cold_ice_complementarity=True,
        normalized_heat_definition='unscaled positive LATENT melt energy divided by each capacity scenario C; K m',
        gain_interpretation='Any later TRAIN gain is empirical output correction, not new physical C or energy-conserving scaled mass.',
        original_v1_failure_preserved=True,field_accuracy_established=False,novelty_claim=False)
    with target.open('x') as stream:stream.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return dict(output=str(target),sha256=sha(target),scenario_methods=len(results),camera_outcomes_read=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['freeze','predict']);p.add_argument('path',type=Path)
    p.add_argument('--out',type=Path);p.add_argument('--commit');a=p.parse_args()
    try:
        if a.mode=='freeze':
            if not a.out:p.error('--out NEW required')
            prepare(a.path,a.out);print('version2 freeze written; no camera or new field predictions read')
        else:print(json.dumps(run(a.path,a.commit or '')))
    except (ValueError,TypeError,KeyError) as error:p.error(str(error))
