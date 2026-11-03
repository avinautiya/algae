"""Additive fixed-alpha phase-model mesh/time/deep-boundary sensitivity.

No camera file, camera outcome or calibration artifact is imported. An optional
explicit empirical output multiplier can be supplied later by the parent; its
scaled retreat is an algebraic diagnostic, not a newly conserved physical solve.
"""
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
from cold_ice_column import simulate
from piramide_experiment import process
from temperature_model import interpolate

ROOT=Path(__file__).resolve().parent
METHODS=('training_temperature_forward_fit','training_derivative_OLS_NOT_source_NTLS')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def resolved(value):
    p=Path(value);return p if p.is_absolute() else ROOT/p


def retreat(melt_mmwe,*,debris_depth_m,ice_probe_gap_m,ice_cell_width_m,ice_depth_m,multiplier=None):
    if type(melt_mmwe) not in (int,float) or not math.isfinite(melt_mmwe) or melt_mmwe<0:raise ValueError('nonnegative finite melt required')
    for value in (debris_depth_m,ice_probe_gap_m,ice_cell_width_m,ice_depth_m):
        if type(value) not in (int,float) or not math.isfinite(value) or value<=0:raise ValueError('positive finite length required')
    displacement=melt_mmwe/917.
    result=dict(raw_ice_equivalent_retreat_m=displacement,
                raw_retreat_over_debris_domain_depth=displacement/debris_depth_m,
                raw_retreat_over_initial_DT5_ice_gap=displacement/ice_probe_gap_m,
                raw_retreat_over_ice_cell_width=displacement/ice_cell_width_m,
                raw_retreat_over_ice_domain_depth=displacement/ice_depth_m,
                retreat_was_not_applied_to_solver_geometry=True)
    if multiplier is not None:
        if type(multiplier) not in (int,float) or not math.isfinite(multiplier) or multiplier<=0:raise ValueError('positive finite empirical multiplier required')
        result['empirically_scaled_diagnostic']=dict(multiplier=multiplier,
            ice_equivalent_retreat_m=displacement*multiplier,
            retreat_over_debris_domain_depth=displacement*multiplier/debris_depth_m,
            retreat_over_initial_DT5_ice_gap=displacement*multiplier/ice_probe_gap_m,
            not_a_conserved_physical_solve=True)
    else:result['empirical_multiplier_not_supplied']=True
    return result


def compare(train,evaluation,original):
    if train<=0 or evaluation<0:raise ValueError('positive train/nonnegative evaluation melt required')
    ratio=evaluation/train
    return dict(raw_training_melt_mmwe=train,raw_evaluation_melt_mmwe=evaluation,
                evaluation_training_raw_melt_ratio=ratio,
                relative_change_raw_training_melt=train/original['raw_training_melt_mmwe']-1,
                relative_change_raw_evaluation_melt=evaluation/original['raw_evaluation_melt_mmwe']-1,
                fixed_original_empirical_gain_prediction_multiplier=evaluation/original['raw_evaluation_melt_mmwe'],
                same_train_target_recalibrated_prediction_multiplier=ratio/(original['raw_evaluation_melt_mmwe']/original['raw_training_melt_mmwe']),
                relative_change_raw_melt_ratio=ratio/(original['raw_evaluation_melt_mmwe']/original['raw_training_melt_mmwe'])-1)


def run(normalized_path,freeze_path,output,empirical_component_multiplier=None):
    if output.exists():raise ValueError('refusing to overwrite numerical sensitivity')
    p=json.loads(normalized_path.read_text());f=json.loads(freeze_path.read_text());c=f['config']
    if p['freeze_sha256']!=sha(freeze_path):raise ValueError('phase artifact does not match freeze')
    for key,expected in f['source_hashes'].items():
        if sha(ROOT/key)!=expected:raise ValueError('frozen phase source changed')
    for key,expected in f['input_hashes'].items():
        if sha(resolved(c[key]))!=expected:raise ValueError('frozen provenance input changed')
    v1=json.loads(resolved(c['v1_freeze']).read_text())
    for key,expected in f['v1_process_input_hashes'].items():
        if sha(resolved(v1['config'][key]))!=expected:raise ValueError('frozen process input changed')
    times,data=process(resolved(v1['config']['process_csv']),resolved(v1['config']['process_manifest']))
    start,split,end=[dt.datetime.fromisoformat(p[k]) for k in ('start_local','split_local','end_local')]
    chosen=[i for i,t in enumerate(times) if start<=t<=end];labels=[times[i] for i in chosen];values=[data[i] for i in chosen]
    if labels[0]!=start or labels[-1]!=end or split not in labels:raise ValueError('exact source-clock endpoint support required')
    cut=labels.index(split);results=[]
    cases=[(2*c['debris_cells'],2*c['ice_cells'],c['dt_seconds']/2,c['ice_depth_m']),
           (2*c['debris_cells'],2*c['ice_cells'],c['dt_seconds']/4,c['ice_depth_m']),
           (2*c['debris_cells'],4*c['ice_cells'],c['dt_seconds']/2,2*c['ice_depth_m'])]
    for method in METHODS:
        rows=[r for r in p['results'] if r['method']==method and r['geometry']['id']==p['primary_geometry']
              and r['hold']==p['primary_hold'] and r['capacity_id']==p['primary_capacity_id']]
        if len(rows)!=1 or rows[0]['status']!='predicted':raise ValueError('one successfully predicted primary-capacity method required')
        base=rows[0];geom=base['geometry'];top=geom['sensor_depths_m'][0]
        depths=[z-top for z in geom['sensor_depths_m']];depth=geom['ice_depth_m']-top;gap=geom['ice_depth_m']-geom['sensor_depths_m'][-1]
        original_total=base['raw_training_melt_mmwe']+base['raw_evaluation_melt_mmwe']
        result=dict(method=method,geometry=geom,hold=base['hold'],capacity_id=base['capacity_id'],
            capacity_J_m3K=base['physicalcapacity_J_m3K'],alpha_m2_s=base['alpha_m2_s'],
            original=dict(debris_cells=c['debris_cells'],ice_cells=c['ice_cells'],dt_seconds=c['dt_seconds'],ice_depth_m=c['ice_depth_m'],
                          raw_training_melt_mmwe=base['raw_training_melt_mmwe'],raw_evaluation_melt_mmwe=base['raw_evaluation_melt_mmwe'],
                          energy_ledger=base['physical_energy_ledger'],
                          geometry_diagnostic=retreat(original_total,debris_depth_m=depth,ice_probe_gap_m=gap,
                              ice_cell_width_m=c['ice_depth_m']/c['ice_cells'],ice_depth_m=c['ice_depth_m'],
                              multiplier=empirical_component_multiplier)),cases=[])
        for nd,ni,step,ice_depth in cases:
            initial=interpolate(depths+[depth],values[0]+[0.],[(i+.5)*depth/nd for i in range(nd)])
            upper=[row[0] for row in (values[:-1] if base['hold']=='left' else values[1:])]
            solved=simulate(debris_conductivity_W_mK=base['debris_conductivity_W_mK'],
                debris_heat_capacity_J_m3K=base['physicalcapacity_J_m3K'],debris_depth_m=depth,
                initial_debris_temperature_C=initial,initial_ice_temperature_C=[c['cold_ice_initial_temperature_C']]*ni,
                surface_temperature_C=upper,interval_seconds=1800.,dt_seconds=step,ice_depth_m=ice_depth,
                ice_conductivity_W_mK=c['ice_conductivity_W_mK'],ice_heat_capacity_J_m3K=c['ice_heat_capacity_J_m3K'],
                latent_heat_J_kg=c['latent_heat_J_kg'])
            intervals=solved['intervals'];train=math.fsum(x['melt_mmwe'] for x in intervals[:cut]);evaluation=math.fsum(x['melt_mmwe'] for x in intervals[cut:])
            row=compare(train,evaluation,base)
            history=solved['history'];row.update(debris_cells=nd,ice_cells=ni,dt_seconds=step,ice_depth_m=ice_depth,
                relative_energy_residual=solved['relative_energy_residual'],energy_residual_J_m2=solved['energy_residual_J_m2'],
                maximum_step_energy_residual_J_m2=solved['max_step_energy_residual_J_m2'],
                complementarity_max_product_C_kg_m2=solved['complementarity_max_product_C_kg_m2'],
                minimum_interface_temperature_C=min(x['interface_temperature_C'] for x in history),
                minimum_melt_substep_kg_m2=min(x['melt_kg_m2'] for x in history),
                deep_heat_out_J_m2=solved['deep_heat_out_J_m2'],surface_heat_J_m2=solved['surface_heat_J_m2'],
                deep_out_heat_over_surface_heat=(solved['deep_heat_out_J_m2']/solved['surface_heat_J_m2'] if solved['surface_heat_J_m2'] else None),
                final_ice_temperature_min_C=min(solved['final_ice_temperature_C']),
                final_ice_temperature_max_C=max(solved['final_ice_temperature_C']),
                geometry_diagnostic=retreat(train+evaluation,debris_depth_m=depth,ice_probe_gap_m=gap,
                    ice_cell_width_m=ice_depth/ni,ice_depth_m=ice_depth,multiplier=empirical_component_multiplier))
            result['cases'].append(row)
        a,b,z=result['cases'];result['paired_differences']=dict(
            time_dt150_to_dt75_relative_melt_ratio_change=b['evaluation_training_raw_melt_ratio']/a['evaluation_training_raw_melt_ratio']-1,
            ice_depth2_to4m_same_cell_width_relative_melt_ratio_change=z['evaluation_training_raw_melt_ratio']/a['evaluation_training_raw_melt_ratio']-1)
        results.append(result)
    report=dict(normalized_predictions_sha256=sha(normalized_path),phase_freeze_sha256=sha(freeze_path),
        refinement_source_sha256=sha(Path(__file__)),results=results,camera_data_imported=False,alpha_refitted=False,
        initial_ice_temperature_changed=False,original_primary_prediction_replaced=False,
        empirical_component_multiplier=empirical_component_multiplier,
        scope='conditional fixed-parameter mesh/time/deep-boundary sensitivities, not certified error bounds or statistical confidence intervals')
    with output.open('x') as stream:stream.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('normalized',type=Path);p.add_argument('freeze',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--empirical-component-multiplier',type=float);a=p.parse_args()
    try:
        result=run(a.normalized,a.freeze,a.output,a.empirical_component_multiplier)
        print(json.dumps(dict(output=str(a.output),methods=len(result['results']),camera_data_imported=False)))
    except (ValueError,TypeError,KeyError) as error:p.error(str(error))
