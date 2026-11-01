"""Strict self-reported evidence checklist, never an independent validation."""
import hashlib
import json
import math

GATES=('hydraulic_partition','water_conservation','thermal_latent_conservation',
       'spatial_convergence','temporal_convergence','event_timing',
       'closure_sensitivity','boundary_applicability')


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def sha(value):
    return isinstance(value,str) and len(value)==64 and all(c in '0123456789abcdef' for c in value)


def readiness(record):
    """Require explicit gates and provenance; a true return is not scientific approval.

    This checks declared records, not their truth. External validation is an
    additional gate for an environmental-performance claim. Publication quality
    and novelty are deliberately not machine-approved.
    """
    if not isinstance(record,dict):raise ValueError('study record must be an object')
    # Reject NaN/Infinity even in fields not otherwise inspected.
    fingerprint=digest(record)
    problems=[]
    scope=record.get('claim_scope')
    if scope not in ('conditional_computational','environmental_prediction'):
        problems.append('claim_scope must be conditional_computational or environmental_prediction')
    if not isinstance(record.get('study_id'),str) or not record.get('study_id'):
        problems.append('missing study_id')
    freeze=record.get('freeze',{})
    if not isinstance(freeze,dict):freeze={}
    if not sha(freeze.get('configuration_sha256')):problems.append('missing valid frozen configuration hash')
    if freeze.get('predictions_frozen_before_test') is not True:
        problems.append('predictions were not declared frozen before test')
    commit=freeze.get('commit')
    if not isinstance(commit,str) or len(commit)!=40 or any(c not in '0123456789abcdef' for c in commit):
        problems.append('missing pinned full commit')
    gates=record.get('gates',{})
    if not isinstance(gates,dict):gates={}
    for name in GATES:
        gate=gates.get(name,{})
        if not isinstance(gate,dict) or gate.get('status')!='pass' or not sha(gate.get('record_sha256')):
            problems.append(name+': absent, not passed, or missing evidence hash')
    campaign=record.get('campaign',{})
    if not isinstance(campaign,dict):campaign={}
    expected=campaign.get('expected_case_ids');done=campaign.get('completed_case_ids')
    def ids(x):return isinstance(x,list) and all(isinstance(v,str) and v for v in x) and len(set(x))==len(x) and bool(x)
    if not ids(expected) or not ids(done) or set(expected)!=set(done):
        problems.append('complete unique expected/completed case sets must agree')
    if campaign.get('unresolved_failure_ids')!=[]:problems.append('failure ledger missing or unresolved failures remain')
    if scope=='environmental_prediction':
        val=record.get('external_validation',{})
        if not isinstance(val,dict) or val.get('independent') is not True or val.get('status')!='pass' or not sha(val.get('record_sha256')):
            problems.append('independent external validation missing, not passed, or undocumented')
    return dict(checklist_complete=not problems,problems=problems,claim_scope=scope,
                record_sha256=fingerprint,
                qualification='Checks self-reported metadata only; not independent verification, priority, publication or real-glacier approval')
