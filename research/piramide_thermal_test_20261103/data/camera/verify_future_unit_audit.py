"""Check already role-released endpoint audit only; never reads native MAT files."""
from pathlib import Path
import hashlib,json,math
P=Path(__file__).parent
source=P/'FUTURE_BATCH_ENDPOINT_SOURCE_AUDIT.json';r=json.loads(source.read_text())
expected=[('2022-12-30T00:00:00','2023-01-13T03:00:00'),('2023-01-13T03:00:00','2023-01-27T03:00:00'),('2023-01-27T03:00:00','2023-02-10T03:00:00')]
assert [(x['start'],x['end']) for x in r['intervals']]==expected
assert len(r['endpoints'])==4 and not any('2023-02-24' in x['source_local_label'] for x in r['endpoints'])
rows=[]
for x in r['intervals']:
    native=x['independent_ice_conversion_mmwe'];provided=x['provided_ablation_change_source_units_unverified']*1000
    assert math.isclose(native,x['raw_distance_change_m_ice']*917,rel_tol=1e-13,abs_tol=1e-12)
    assert math.isclose(native,provided,rel_tol=1e-13,abs_tol=1e-12)
    assert x['operator_changed'] is False
    rows.append({'start':x['start'],'end':x['end'],'raw_distance_delta_times917_mmwe':native,
                 'provided_ablation_delta_times1000_mmwe':provided,'absolute_conversion_difference_mmwe':abs(native-provided),
                 'operator_changed':False,'negative_endpoint_difference':x['negative_distance_change']})
out={'role_release_audit_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
     'source_member_sha256':r['source_sha256'],'interval_unit_checks':rows,
     'excluded_Feb24_value_accessed':False,'passed':True,
     'claim':'Unit/operator check, not proof of camera precision, reference continuity, independent glacier replication or predictive superiority',
     'camera_endpoint_precision':'Unknown; no numerical delta or confidence interval justified',
     'manual_corrections_and_resets':'Only four authorized endpoints audited; interior changes unresolved; no operator switched'}
(P/'FUTURE_BATCH_UNIT_VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
print('Three released future interval conversions verified; excluded endpoint not accessed.')
