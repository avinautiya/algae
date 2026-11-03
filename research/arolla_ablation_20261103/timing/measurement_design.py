"""Retrospective resolution diagnostic, not a prospective power calculation."""
import itertools,json,math
from pathlib import Path

def reading_resolution(fa,fb,train,test):
    """Critical per-reading bound for sign of shared-reading squared-error gap.

    Fixed driver ratios and positive true training increment are required.
    Strict inequality delta < reported threshold is needed; equality permits ties.
    """
    if not all(math.isfinite(x) for x in (fa,fb,train,test)) or min(fa,fb)<0 or train<=0 or test<0:
        raise ValueError('invalid ratios/increments')
    if fa==fb: return {'status':'identical_predictions','critical_bound_mmwe':0.0}
    midpoint=(fa+fb)*train-2*test
    critical=min(abs(midpoint)/(2*(fa+fb)+4),train/2,test/2)
    gap=(fa-fb)*midpoint
    return {'status':'strict_nominal_ranking' if gap else 'nominal_tie',
            'nominal_better':'a' if gap<0 else 'b' if gap>0 else None,
            'critical_bound_mmwe':critical,
            'critical_per_reading_ice_bound_m_at_rho900':critical/900,
            'strict_inequality_required':True}

def main():
    root=Path(__file__).resolve().parents[1]
    data=json.loads((root/'model/results/frozen_ablation_comparison/RESULTS.json').read_text())
    groups={}
    for r in data['rows']:
        if r['stake']!='H1D' or 'test_train_dose_ratio' not in r: continue
        key=(r['mean_label_end_offset_h'],tuple(r['visit_hours']))
        groups.setdefault(key,{})[r['model']]=r
    pairs=[]
    for a,b in itertools.combinations(sorted(next(iter(groups.values()))),2):
        entries=[]
        for clock,g in groups.items():
            ra,rb=g[a],g[b]
            v=reading_resolution(ra['test_train_dose_ratio'],rb['test_train_dose_ratio'],ra['calibration_observed_mmwe'],ra['heldout_observed_mmwe'])
            entries.append(dict(clock=[clock[0],list(clock[1])],nominal=ra['nominal'],**v))
        winners={x['nominal_better'] for x in entries}
        pairs.append(dict(a=a,b=b,clock_only_ranking_consistent=len(winners)==1,
                          nominal=next(x for x in entries if x['nominal']),
                          smallest_conditional_threshold_mmwe=min(x['critical_bound_mmwe'] for x in entries)))
    out=dict(pairs=pairs,retrospective_design_diagnostic=True,new_field_accuracy_evidence=False,
             limitation='Thresholds concern fixed driver ratios and conditional bounded reading errors. They do not resolve model discrepancy, snow, unknown error semantics or uncertain visit time.')
    (root/'timing/MEASUREMENT_DESIGN.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
