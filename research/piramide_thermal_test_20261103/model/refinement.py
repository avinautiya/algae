"""Additive fixed-parameter refinement audit; never imports camera outcomes.

This diagnostic does not amend the frozen primary prediction or selection rule.
Both a fixed melt-gain interpretation and reapplying the same TRAIN-only gain
rule can be evaluated without knowing a camera target: their multipliers depend
only on normalized thermal integrals.
"""
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
from piramide_experiment import process
from temperature_model import predict

CASES=((32,150.),(64,150.),(64,75.))
METHODS=('training_temperature_forward_fit','training_derivative_OLS_NOT_source_NTLS')
ROOT=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metrics(heat,split_index,baseline):
    train=math.fsum(heat[:split_index]);evaluation=math.fsum(heat[split_index:])
    if train<=0 or evaluation<=0:raise ValueError('positive train/evaluation heat required for ratio diagnostic')
    ratio=evaluation/train
    return dict(training_normalized_basal_heat_K_m=train,
                evaluation_normalized_basal_heat_K_m=evaluation,
                evaluation_training_heat_ratio=ratio,
                fixed_original_melt_gain_prediction_multiplier=evaluation/baseline['evaluation_normalized_basal_heat_K_m'],
                same_train_target_recalibrated_gain_prediction_multiplier=ratio/baseline['evaluation_training_heat_ratio'],
                relative_change_train_heat=train/baseline['training_normalized_basal_heat_K_m']-1,
                relative_change_evaluation_heat=evaluation/baseline['evaluation_normalized_basal_heat_K_m']-1,
                relative_change_heat_ratio=ratio/baseline['evaluation_training_heat_ratio']-1,
                prediction_multiplier_units='dimensionless; no camera target needed')


def run(normalized_path,freeze_path,output):
    if output.exists():raise ValueError('refusing to overwrite refinement diagnostic')
    p=json.loads(normalized_path.read_text());f=json.loads(freeze_path.read_text());c=f['config']
    if p['freeze_sha256']!=sha(freeze_path):raise ValueError('normalized predictions do not match freeze')
    # All initial frozen inputs/sources must remain intact. Resolve the source
    # input paths relative to model/, as in the published campaign commands.
    for name,expected in f['source_hashes'].items():
        if sha(ROOT/name)!=expected:raise ValueError('frozen source changed: '+name)
    for name,expected in f['input_hashes'].items():
        path=Path(c[name]);path=path if path.is_absolute() else ROOT/path
        if sha(path)!=expected:raise ValueError('frozen input changed: '+name)
    def local(name):
        path=Path(c[name]);return path if path.is_absolute() else ROOT/path
    times,values=process(local('process_csv'),local('process_manifest'))
    start,split,end=[dt.datetime.fromisoformat(p[name]) for name in ('start_local','split_local','end_local')]
    indices=[i for i,t in enumerate(times) if start<=t<=end]
    labels=[times[i] for i in indices];v=[values[i] for i in indices]
    if labels[0]!=start or labels[-1]!=end or split not in labels:raise ValueError('exact native endpoint support required')
    split_index=labels.index(split);results=[]
    for method in METHODS:
        matches=[r for r in p['results'] if r['geometry']['id']==p['primary_geometry'] and r['hold']==p['primary_hold'] and r['method']==method]
        if len(matches)!=1:raise ValueError('exactly one frozen primary method record required')
        original=matches[0]
        if original['status']!='predicted':
            results.append(dict(method=method,status='frozen_method_refused',reason=original.get('reason')));continue
        alpha=original['alpha_m2_s'];geometry=original['geometry'];top=geometry['sensor_depths_m'][0]
        sensor_depths=[z-top for z in geometry['sensor_depths_m']];depth=geometry['ice_depth_m']-top
        rows=[]
        for cells,step in CASES:
            r=predict(alpha,domain_depth_m=depth,initial_depths_m=sensor_depths+[depth],
                initial_temperature_C=v[0]+[0.],upper_temperature_C=[x[0] for x in v],
                lower_temperature_C=[0.]*len(v),interval_seconds=1800.,
                sample_depths_m=sensor_depths[1:],cells=cells,dt_seconds=step,hold=original['hold'])
            heat=r['interval_basal_normalized_heat_K_m'];q=metrics(heat,split_index,original)
            scale=max(abs(q['training_normalized_basal_heat_K_m'])+abs(q['evaluation_normalized_basal_heat_K_m']),1.)
            q.update(cells=cells,dt_seconds=step,alpha_m2_s=alpha,
                normalized_energy_residual_K_m=r['normalized_energy_residual_K_m'],
                relative_energy_residual_to_basal_integral=r['normalized_energy_residual_K_m']/scale,
                max_absolute_step_normalized_energy_residual_K_m=r['maximum_step_normalized_energy_residual_K_m'],
                negative_basal_flux_substeps=r['negative_basal_flux_substeps'],
                negative_basal_heat_intervals=sum(x<0 for x in heat),
                no_refreezing_melt_conversion_eligible=r['negative_basal_flux_substeps']==0 and all(x>=0 for x in heat))
            rows.append(q)
        paired=dict(spatial_N32_to_N64_at_dt150_heat_ratio_relative_change=rows[1]['evaluation_training_heat_ratio']/rows[0]['evaluation_training_heat_ratio']-1,
                    temporal_dt150_to_dt75_at_N64_heat_ratio_relative_change=rows[2]['evaluation_training_heat_ratio']/rows[1]['evaluation_training_heat_ratio']-1)
        results.append(dict(method=method,status='refined',geometry=geometry,hold=original['hold'],selected_alpha_fixed=True,
            original_cells=c['cells'],original_dt_seconds=c['dt_seconds'],original_alpha_m2_s=alpha,
            original_training_normalized_basal_heat_K_m=original['training_normalized_basal_heat_K_m'],
            original_evaluation_normalized_basal_heat_K_m=original['evaluation_normalized_basal_heat_K_m'],
            original_evaluation_training_heat_ratio=original['evaluation_training_heat_ratio'],cases=rows,paired_differences=paired))
    report=dict(normalized_predictions_sha256=sha(normalized_path),original_freeze_sha256=sha(freeze_path),
        refinement_source_sha256=sha(Path(__file__)),source_clock_start=start.isoformat(),source_clock_split=split.isoformat(),
        source_clock_end=end.isoformat(),results=results,camera_outcomes_loaded=False,alpha_refitted=False,
        original_study_criterion_changed=False,original_primary_prediction_replaced=False,
        scope='fixed-parameter mesh/time sensitivity; not a rigorous numerical error bound or statistical confidence interval')
    with output.open('x') as stream:stream.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('normalized',type=Path);p.add_argument('freeze',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    try:
        report=run(a.normalized,a.freeze,a.output)
        print(json.dumps(dict(output=str(a.output),methods=len(report['results']),camera_outcomes_loaded=False)))
    except (ValueError,KeyError,TypeError) as error:p.error(str(error))
