"""Long-trajectory fixed-parameter numerical sensitivity, no future truth.

Per-period dose changes are directly the prediction changes under the original
unchanged version3 gains. This module does not read those gains or camera values.
"""
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
from cold_ice_column import simulate
from piramide_experiment import process
from temperature_model import interpolate
from batch_phase_experiment import sum_interval
from phase_refinement import retreat

ROOT=Path(__file__).resolve().parent


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def path(value):
    p=Path(value);return p if p.is_absolute() else ROOT/p


def compare_block(refined,base):
    if refined['start']!=base['start'] or refined['end']!=base['end']:raise ValueError('identical interval operators required')
    denominator=base['normalized_basal_heat_K_m']
    if denominator<=0:raise ValueError('positive original dose required for relative sensitivity')
    q=dict(refined)
    q.update(original_normalized_basal_heat_K_m=denominator,
             original_raw_melt_mmwe=base['raw_melt_mmwe'],
             fixed_version3_gain_prediction_multiplier=refined['normalized_basal_heat_K_m']/denominator,
             relative_fixed_version3_gain_prediction_change=refined['normalized_basal_heat_K_m']/denominator-1)
    return q


def run(normalized,freeze,output):
    if output.exists():raise ValueError('refusing to overwrite long-trajectory diagnostic')
    p=json.loads(normalized.read_text());f=json.loads(freeze.read_text());c=f['config'];physical=f['physical_config']
    if p['freeze_sha256']!=sha(freeze):raise ValueError('prediction/freeze mismatch')
    for name,h in f['source_hashes'].items():
        if sha(ROOT/name)!=h:raise ValueError('frozen model source changed')
    for name,h in f['input_hashes'].items():
        if sha(path(c[name]))!=h:raise ValueError('frozen batch provenance changed')
    v1=json.loads(path(physical['v1_freeze']).read_text())
    for name,h in f['source_v3_process_hashes'].items():
        if sha(path(v1['config'][name]))!=h:raise ValueError('frozen process input changed')
    times,rows=process(path(v1['config']['process_csv']),path(v1['config']['process_manifest']))
    start=dt.datetime.fromisoformat(p['continuous_model_initialization_local'])
    end=dt.datetime.fromisoformat(p['continuous_model_end_local'])
    selected=[i for i,t in enumerate(times) if start<=t<=end];labels=[times[i] for i in selected];data=[rows[i] for i in selected]
    if labels[0]!=start or labels[-1]!=end:raise ValueError('exact full native support required')
    primary=c['primary_candidate'];matched=primary.replace('forward__','derivative_OLS__',1)
    cases=[dict(debris_cells=2*physical['debris_cells'],ice_cells=2*physical['ice_cells'],
                dt_seconds=physical['dt_seconds']/2,ice_depth_m=physical['ice_depth_m']),
           dict(debris_cells=2*physical['debris_cells'],ice_cells=4*physical['ice_cells'],
                dt_seconds=physical['dt_seconds']/2,ice_depth_m=2*physical['ice_depth_m'])]
    results=[]
    for ident in [primary,matched]:
        matches=[r for r in p['results'] if r['model_id']==ident]
        if len(matches)!=1 or matches[0]['status']!='predicted':raise ValueError('exact successfully predicted primary/matched model required')
        old=matches[0];geom=old['geometry'];top=geom['sensor_depths_m'][0]
        zz=[z-top for z in geom['sensor_depths_m']];depth=geom['ice_depth_m']-top;gap=geom['ice_depth_m']-geom['sensor_depths_m'][-1]
        upper=[row[0] for row in (data[:-1] if old['hold']=='left' else data[1:])]
        refined=[]
        for case in cases:
            n=case['debris_cells'];C=old['physicalcapacity_J_m3K'];L=physical['latent_heat_J_kg']
            Td=interpolate(zz+[depth],data[0]+[0.],[(j+.5)*depth/n for j in range(n)])
            result=simulate(debris_conductivity_W_mK=old['debris_conductivity_W_mK'],debris_heat_capacity_J_m3K=C,
                debris_depth_m=depth,initial_debris_temperature_C=Td,
                initial_ice_temperature_C=[physical['cold_ice_initial_temperature_C']]*case['ice_cells'],
                surface_temperature_C=upper,interval_seconds=1800.,dt_seconds=case['dt_seconds'],ice_depth_m=case['ice_depth_m'],
                ice_conductivity_W_mK=physical['ice_conductivity_W_mK'],ice_heat_capacity_J_m3K=physical['ice_heat_capacity_J_m3K'],
                latent_heat_J_kg=L)
            blocks=[]
            for base in old['future_blocks']:
                block=sum_interval(result['intervals'],labels,base['start'],base['end'],C,L)
                blocks.append(compare_block(block,base))
            alpha_ice=physical['ice_conductivity_W_mK']/physical['ice_heat_capacity_J_m3K']
            skin=math.sqrt(alpha_ice*86400/math.pi)
            velocity=(result['melt_mmwe']/917)/(end-start).total_seconds()
            refined.append(dict(numerics=case,future_blocks=blocks,
                cumulative_raw_melt_mmwe=result['melt_mmwe'],
                energy_ledger=dict(initial_sensible_storage_J_m2=result['initial_sensible_storage_J_m2'],
                    final_sensible_storage_J_m2=result['final_sensible_storage_J_m2'],surface_heat_J_m2=result['surface_heat_J_m2'],
                    deep_heat_out_J_m2=result['deep_heat_out_J_m2'],latent_heat_J_m2=result['latent_heat_J_m2'],
                    energy_residual_J_m2=result['energy_residual_J_m2'],relative_energy_residual=result['relative_energy_residual'],
                    max_step_energy_residual_J_m2=result['max_step_energy_residual_J_m2']),
                complementarity_max_product_C_kg_m2=result['complementarity_max_product_C_kg_m2'],
                negative_melt_substeps=sum(row['melt_kg_m2']<0 for row in result['history']),
                minimum_interface_temperature_C=min(row['interface_temperature_C'] for row in result['history']),
                geometry_diagnostic=retreat(result['melt_mmwe'],debris_depth_m=depth,ice_probe_gap_m=gap,
                    ice_cell_width_m=case['ice_depth_m']/case['ice_cells'],ice_depth_m=case['ice_depth_m']),
                omitted_ice_advection_diagnostic=dict(raw_mean_retreat_velocity_m_s=velocity,
                    ice_diffusivity_m2_s=alpha_ice,diurnal_diffusive_skin_depth_m=skin,
                    mean_retreat_peclet_on_diurnal_skin=velocity*skin/alpha_ice,
                    scope='raw mean scale diagnostic only; not a bound on transient or empirically scaled advection error')))
        results.append(dict(model_id=ident,alpha_m2_s=old['alpha_m2_s'],capacity_J_m3K=old['physicalcapacity_J_m3K'],
            original_numerics={k:physical[k] for k in ('debris_cells','ice_cells','dt_seconds','ice_depth_m')},
            original_future_blocks=old['future_blocks'],original_energy_ledger=old['physical_energy_ledger'],
            cases=refined))
    report=dict(normalized_predictions_sha256=sha(normalized),batch_freeze_sha256=sha(freeze),
        numerical_diagnostic_source_sha256=sha(Path(__file__)),continuous_start_local=start.isoformat(),continuous_end_local=end.isoformat(),
        results=results,camera_data_imported=False,gain_values_imported=False,alpha_refitted=False,source_models_changed=False,
        source_study_criterion_changed=False,
        scope='fixed-parameter long-trajectory numerical/deep-boundary sensitivities; not rigorous error bounds, field errors or statistical confidence intervals')
    with output.open('x') as stream:stream.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('normalized',type=Path);p.add_argument('freeze',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    try:
        q=run(a.normalized,a.freeze,a.output);print(json.dumps(dict(output=str(a.output),models=len(q['results']),future_truth_imported=False)))
    except (ValueError,KeyError,TypeError) as error:p.error(str(error))
