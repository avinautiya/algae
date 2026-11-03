"""Deterministic source-missingness amendment; never reads outcome values."""
import copy,datetime as dt,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    directory=ROOT/'model/frozen_piramide_v3'
    if directory.exists():raise ValueError('refuse overwriting coverage amendment')
    source=ROOT/'model/frozen_piramide_v2/FREEZE.json';f=copy.deepcopy(json.loads(source.read_text()))
    audit_path=ROOT/'data/camera/MISSINGNESS_ONLY_ENDPOINT_AUDIT.json';a=json.loads(audit_path.read_text())
    if a['new_target_values_exposed'] is not False or a['prediction_error_computed'] is not False:raise ValueError('coverage-only amendment required')
    old=dt.datetime.fromisoformat(f['config']['end'])
    if old!=dt.datetime.fromisoformat(a['original_end']):raise ValueError('original endpoint mismatch')
    finite=[dt.datetime.fromisoformat(x['source_local_label']) for x in a['nearby_timestamps_finiteness_only'] if x['distance_isfinite'] and abs((dt.datetime.fromisoformat(x['source_local_label'])-old).total_seconds())<=43200]
    if not finite:raise ValueError('no supported endpoint within frozen twelve-hour radius')
    selected=min(finite,key=lambda x:(abs((x-old).total_seconds()),x))
    if selected!=dt.datetime.fromisoformat(a['selected_end']):raise ValueError('nearest-finite rule mismatch')
    end=selected.isoformat();f['config']['end']=end
    for x in f['config']['camera_evaluation']['intervals']:
        if x['role']=='evaluation':x['end']=end
    f['camera_evaluation']=copy.deepcopy(f['config']['camera_evaluation'])
    spec=json.loads((ROOT/'analysis/ANALYSIS_COLD_SPEC.json').read_text())
    spec['primary_evaluation'][1]=end
    spec['limitations'].append('Transparent missingness-only amendment: original endpoint is native NaN; nearest finite endpoint within12h chosen without reading its value or calculating any evaluation error.')
    out_spec=ROOT/'analysis/ANALYSIS_COVERAGE_SPEC.json'
    with out_spec.open('x') as stream:json.dump(spec,stream,indent=2);stream.write('\n')
    f['analysis_spec_sha256']=sha(out_spec)
    f['source_hashes']['prepare_coverage_amendment.py']=sha(Path(__file__))
    f['coverage_amendment']=dict(original_freeze_sha256=sha(source),missingness_audit_sha256=sha(audit_path),
        selected_end=end,original_end=old.isoformat(),rule=a['selection_rule'],
        training_outcomes_previously_seen=True,original_evaluation_value_is_missing=True,
        replacement_evaluation_value_seen=False,evaluation_error_computed=False,
        model_family_parameters_and_primary_unchanged=True)
    f['camera_outcomes_read']=True # TRAIN only + native missing flag; disclose history.
    f['new_field_model_computed']=False
    directory.mkdir()
    (directory/'FREEZE.json').write_text(json.dumps(f,indent=2)+'\n')
    print(json.dumps(dict(selected_end=end,freeze_sha256=sha(directory/'FREEZE.json'),replacement_value_read=False)))
if __name__=='__main__':main()
