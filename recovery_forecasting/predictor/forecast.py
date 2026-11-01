"""Validated JSON CLI for an idealized, static first-contact forecast.

This is a conditional model feasibility calculation, not glacier validation,
not a prescribed-hydrograph forecast, and not sustained ordinary-flow recovery.
"""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import warnings
from scipy.integrate import IntegrationWarning

try:
    from .vendor.recovery_budget import static_bound
except ImportError:
    from vendor.recovery_budget import static_bound

VERSION = '1.0.0'
KERNEL_SHA256 = '43c3a459cd70722c754d9c63aa78f2f80950a7733c096710180546b35289cc78'
PARAMETERS = ('area_m2', 'coefficient', 'exponent', 'bypass_m', 'toe_m',
              'initial_stage_m', 'peak_m3_s')
POSITIVE = {'area_m2', 'coefficient', 'exponent', 'peak_m3_s'}
SCOPE = ('Static reservoir and bypass; freely selectable bounded inflow; first '
         'zero-depth contact only. Not a prescribed flood, positive contact '
         'flux, sustained recovery, cliff lifetime, or field validation.')


class InputError(ValueError):
    pass


def _keys(value, required, optional=(), label='object'):
    if not isinstance(value, dict):raise InputError(f'{label} must be an object')
    missing=set(required)-value.keys();extra=value.keys()-set(required)-set(optional)
    if missing:raise InputError(f'{label}: missing {sorted(missing)}')
    if extra:raise InputError(f'{label}: unknown fields {sorted(extra)}')


def _number(value, label, *, positive=False, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value,(int,float)):
        raise InputError(f'{label} must be a finite number, not a bool/string')
    try:valid=math.isfinite(value)
    except (OverflowError,TypeError):valid=False
    if not valid or (positive and value<=0) or (nonnegative and value<0):
        raise InputError(f'{label} violates its finite/physical range')
    return float(value)


def validate(document):
    _keys(document,('schema_version','evidence_label','parameters','budgets'),
          ('parameter_bounds','notes'),label='document')
    if type(document['schema_version']) is not int or document['schema_version']!=1:
        raise InputError('schema_version must be integer 1')
    if not isinstance(document['evidence_label'],str) or not document['evidence_label'].strip():
        raise InputError('nonempty evidence_label required')
    if 'notes' in document and not isinstance(document['notes'],str):
        raise InputError('notes must be a string')
    _keys(document['parameters'],PARAMETERS,label='parameters')
    params={key:_number(document['parameters'][key],key,positive=key in POSITIVE) for key in PARAMETERS}
    _keys(document['budgets'],('available_volume_m3','available_time_s'),label='budgets')
    budgets={key:_number(value,key,nonnegative=True) for key,value in document['budgets'].items()}
    bounds=document.get('parameter_bounds',{})
    _keys(bounds,(),PARAMETERS,label='parameter_bounds')
    clean_bounds={}
    for key,interval in bounds.items():
        if not isinstance(interval,list) or len(interval)!=2:
            raise InputError(f'parameter_bounds.{key} must be [lower, upper]')
        low,high=(_number(v,key,positive=key in POSITIVE) for v in interval)
        if low>high:raise InputError(f'parameter_bounds.{key} is reversed')
        if not low<=params[key]<=high:raise InputError(f'nominal {key} lies outside bounds')
        clean_bounds[key]=(low,high)
    return params,budgets,clean_bounds


def _classify(params,budgets):
    with warnings.catch_warnings():
        warnings.simplefilter('error',IntegrationWarning)
        b=static_bound(**params)
    initial=params['initial_stage_m']>=params['toe_m']
    if b.reachable_at_peak and (not math.isfinite(b.minimum_time_s) or not math.isfinite(b.minimum_volume_m3)):
        raise ArithmeticError('quadrature returned non-finite reachable prediction')
    if b.minimum_time_s<0 or b.minimum_volume_m3<0:
        raise ArithmeticError('quadrature returned negative prediction')
    failures=[]
    if not initial:
        if not b.reachable_at_peak:failures.append('peak_at_or_below_static_activation_discharge')
        else:
            if budgets['available_time_s']<b.minimum_time_s:failures.append('insufficient_deadline')
            if budgets['available_volume_m3']<b.minimum_volume_m3:failures.append('insufficient_total_inflow_volume')
    classification='initially_in_contact' if initial else ('excluded_static_model' if failures else 'feasible_static_model')
    finite=lambda x:float(x) if math.isfinite(x) else None
    critical=b.critical_peak_m3_s
    return dict(classification=classification,first_contact_feasible=not failures,
                failure_reasons=failures,critical_peak_m3_s=finite(critical),
                minimum_time_s=finite(b.minimum_time_s),minimum_volume_m3=finite(b.minimum_volume_m3),
                time_margin_s=finite(budgets['available_time_s']-b.minimum_time_s),
                volume_margin_m3=finite(budgets['available_volume_m3']-b.minimum_volume_m3),
                peak_margin_m3_s=finite(params['peak_m3_s']-critical),
                numerical_note='No safety margin is inferred. Equality and near-threshold decisions require numerical and parameter error assessment.')


def forecast(document):
    params,budgets,bounds=validate(document)
    nominal=_classify(params,budgets)
    ordered=sorted(bounds)
    corners=[]
    choices=[sorted(set(bounds[key])) for key in ordered]
    for values in itertools.product(*choices):
        if not ordered:break
        p=dict(params);p.update(zip(ordered,values))
        corners.append(dict(parameters=p,prediction=_classify(p,budgets)))
    summary={}
    for key in ('minimum_time_s','minimum_volume_m3','critical_peak_m3_s'):
        vals=[r['prediction'][key] for r in corners]
        finite=[v for v in vals if v is not None]
        summary[key]=dict(minimum_finite_sample=min(finite) if finite else None,
                          maximum_finite_sample=max(finite) if finite else None,
                          nonfinite_or_unreachable_count=sum(v is None for v in vals))
    canonical=json.dumps(document,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    kernel=Path(__file__).parent/'vendor'/'recovery_budget.py'
    actual_kernel_sha=hashlib.sha256(kernel.read_bytes()).hexdigest()
    if actual_kernel_sha!=KERNEL_SHA256:raise InputError('vendored kernel integrity differs from pinned source')
    return dict(schema_version=1,tool='ice-contact-first-contact-forecast',tool_version=VERSION,
                evidence_label=document['evidence_label'],scope=SCOPE,
                provenance=dict(input_canonical_sha256=hashlib.sha256(canonical).hexdigest(),
                                kernel_sha256=actual_kernel_sha,
                                wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                kernel_source_commit='0b8b7565f869dd5ff456ac309b266257fe71faaa'),
                nominal_parameters=params,budgets=budgets,nominal=nominal,
                corner_sensitivity=dict(count=len(corners),
                    interpretation='Deterministic corner samples only; not a confidence interval, probability, or guaranteed interior extrema.',
                    classifications=sorted(set(r['prediction']['classification'] for r in corners)),
                    ranges=summary,samples=corners),
                assumptions=['All water sources counted in total inflow.',
                             'Toe activation elevation and reservoir geometry stay fixed.',
                             'Static monotone bypass law C*max(H-z_b,0)^p.',
                             'No waves, inertia, uncounted storage release, toe lowering or grain failure.',
                             'Parameters and bounds supplied independently; no calibration performed.',
                             'Coefficient units depend on exponent; varied exponent/coefficient corners require physical justification.'])


def load_document(path):
    def pairs(items):
        out={}
        for key,value in items:
            if key in out:raise InputError(f'duplicate JSON key: {key}')
            out[key]=value
        return out
    def bad_constant(value):raise InputError(f'non-standard JSON constant: {value}')
    return json.loads(Path(path).read_text(),object_pairs_hook=pairs,parse_constant=bad_constant)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--output',type=Path)
    args=parser.parse_args(argv)
    try:
        if args.output and args.output.resolve()==args.input.resolve():
            raise InputError('output must not overwrite the input document')
        result=forecast(load_document(args.input))
        text=json.dumps(result,indent=2,allow_nan=False)+'\n'
        if args.output:
            temporary=None
            try:
                with tempfile.NamedTemporaryFile(mode='w',dir=args.output.parent,prefix='.forecast-',delete=False) as handle:
                    temporary=Path(handle.name);handle.write(text);handle.flush();os.fsync(handle.fileno())
                os.replace(temporary,args.output)
            finally:
                if temporary is not None and temporary.exists():temporary.unlink()
        else:sys.stdout.write(text)
    except (InputError,ValueError,ArithmeticError,OverflowError,OSError,IntegrationWarning) as exc:
        sys.stderr.write(json.dumps(dict(error=type(exc).__name__,message=str(exc)),allow_nan=False)+'\n')
        return 2
    return 0


if __name__=='__main__':raise SystemExit(main())
