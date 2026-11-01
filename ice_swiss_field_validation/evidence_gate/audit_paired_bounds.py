"""Reproduce conditional paired bounds from an existing saved result JSON."""
import argparse
import hashlib
import json
from pathlib import Path
from paired_evidence import compare


def audit(path, scenario):
    data = json.loads(path.read_text())
    result = {}
    for label, use_range in [('nominal', False), ('endpoint_sum_engineering_scenario', True)]:
        rows = []
        for r in data['scenarios'][scenario]['accepted_intervals']:
            if r['role'] != 'test':
                continue
            y = r['observed_lowering_m']
            e = r['summed_epoch_error_scenario_m'] if use_range else 0.
            rows.append(dict(id=r['id'], history=r['station'], target_type='direct_bed_lowering',
                unit='m', observation=[y-e, y+e],
                predictions={k: [v, v] for k, v in r['predictions_m'].items()}))
        result[label] = compare(rows, 'db04_scale_fitted',
                               ['constant_h_fitted', 'gnielinski_smooth_unfitted'])
    result['source_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    result['caveat'] = ('Endpoint-summed source standard deviations are used ONLY as analyst '
        'engineering ranges, not hard truth bounds or confidence intervals. Shared systematic '
        'uncertainty makes independently extremized targets conservative. P3/P5 share one channel event.')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('results', type=Path)
    parser.add_argument('--scenario', default='clock+0h_T+0.00C_D1.30m')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(audit(args.results, args.scenario), indent=2)+'\n')
