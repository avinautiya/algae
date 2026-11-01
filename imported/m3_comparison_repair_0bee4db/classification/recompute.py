"""Reclassify frozen PREDICTIONS; CFD labels below are report-derived only."""
import hashlib
import json
from pathlib import Path
from recovery import classify,history_key


ROOT=Path(__file__).resolve().parent
SNAP=ROOT.parent/'snapshot'
PIN='f3602e35e1595bcdbee4208a354a34944f1ae974'


def main():
    spec=json.loads((SNAP/'frozen/production_spec.json').read_text())
    report=(SNAP/'REPORT_M3_FINAL.md').read_text()
    reported=[]
    for line in report.splitlines():
        if not line.startswith('| E') or '| complete |' not in line:continue
        cols=[s.strip() for s in line.strip('|').split('|')]
        name=cols[0];definition=next(x for x in spec['evaluation_runs'] if x['name']==name)
        loss=None if cols[2]=='—' else float(cols[2])
        reconnection=None if cols[4]=='—' else float(cols[4])
        rem=definition['removal_day']
        if loss is None:label='maintained_contact'
        elif rem is None:label='lost_contact_no_intervention'
        elif reconnection is not None:label='genuine_reconnection'
        else:label='failed_recovery_by_observed_horizon'
        entry=dict(name=name,layout=definition['layout'],removal_day=rem,loss_day_reported=loss,
                   recovery_day_reported=reconnection,classification_report_derived=label,
                   evidence='published report table, NOT recomputed from raw CFD time series')
        if loss is not None:
            geometry=dict(spec['geometry']);geometry.pop('dx',None)
            initial=dict(geometry=geometry,layout=definition['layout'],initial=spec['physics']['initial'])
            physics={k:v for k,v in spec['physics'].items() if k not in ('Q_in','T_in','initial')}
            interventions=[] if rem is None else [dict(time_day=rem,operation='complete_prescribed_removal')]
            entry['pre_loss_history_group']=history_key(initial_state=initial,physics=physics,
                boundary=spec['boundary'],forcing_history=[dict(time_day=0,Q_in=spec['physics']['Q_in'],T_in=spec['physics']['T_in'])],
                interventions=interventions,event_day=loss)
        reported.append(entry)
    predicted=[]
    for run in spec['evaluation_runs']:
        for pred in spec['reduced_predictions'].get(run['name'],[]):
            rec=classify(pred['series_days'],pred['series_f'],actual_removal_day=run['removal_day'],
                         requested_removal_day=run['removal_day'],expected_horizon_day=run['days'])
            predicted.append(dict(run=run['name'],method=pred['method'],classification=rec,
                                  frozen_full_resolution_loss_day=pred['loss'],
                                  frozen_full_resolution_reconnection_day=pred['sustained_reconnection']))
    groups={}
    for r in reported:
        if r.get('pre_loss_history_group'):groups.setdefault(r['pre_loss_history_group'],[]).append(r['name'])
    removal_rows=[r for r in reported if r['removal_day'] is not None]
    counts={k:sum(r['classification_report_derived']==k for r in removal_rows) for k in
            ('maintained_contact','genuine_reconnection','failed_recovery_by_observed_horizon')}
    result=dict(pinned_commit=PIN,
        source_sha256={name:hashlib.sha256((SNAP/name).read_bytes()).hexdigest() for name in ('frozen/production_spec.json','REPORT_M3_FINAL.md')},
        raw_cfd_records_available_in_snapshot=False,
        report_derived_CFD_classes=reported,post_removal_report_class_counts=counts,
        distinct_reported_loss_histories=dict(number_of_reported_loss_cases=sum(len(v) for v in groups.values()),
                                             distinct_physical_input_history_groups=len(groups),groups=groups,
                                             warning='deterministic shared-history accounting, NOT a statistical independence claim'),
        recomputed_frozen_daily_predictions=predicted,
        warnings=['CFD classes are transcribed from report event summaries; raw time-series are required for independent reclassification.',
                  'Stored frozen prediction series are downsampled to daily observations; bracketed recomputation cannot replace full-resolution frozen event times.',
                  'Loss-time matches in E04/E06/E07 are one shared pre-removal strip history, not three independent loss tests.',
                  'The four nonporous removal configurations contain two maintained-contact cases, one genuine recovery and one failed recovery.'])
    (ROOT/'corrected_accounting.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(removal_classes=counts,distinct_loss_history_groups=len(groups),reported_loss_cases=sum(len(v) for v in groups.values())),indent=2))


if __name__=='__main__':main()
