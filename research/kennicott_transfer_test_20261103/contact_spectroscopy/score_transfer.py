import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'kennicott_protocol'), str(ROOT / 'review')]
from outcome_loader import load_after_freeze, sha
from decision_rule import evaluate

if __name__ == '__main__':
    receipt_path = Path(sys.argv[1])
    receipt = json.loads(receipt_path.read_text())
    spec = json.loads((receipt_path.parent / receipt['frozen_spec_path']).read_text())
    for name, digest in spec['source_hashes'].items():
        if sha(ROOT / name) != digest:
            raise ValueError(f'Frozen source changed: {name}')
    predictions = json.loads((receipt_path.parent / receipt['predictions_path']).read_text())
    truth = load_after_freeze(receipt_path, expected_role='evaluation')['rows']
    metrics, refusals = {}, {}
    for name, model in predictions['models'].items():
        if model['status'] != 'predicted':
            refusals[name] = model
            continue
        if [row['site'] for row in model['rows']] != [row['site'] for row in truth]:
            raise ValueError('Site mismatch')
        errors = [p['predicted_mmwe'] - t['observed_net_ablation_mmwe']
                  for p, t in zip(model['rows'], truth)]
        metrics[name] = dict(MAE_mmwe=sum(map(abs, errors)) / 3,
                             absolute_errors_mmwe=list(map(abs, errors)), signed_errors_mmwe=errors)
    decision = evaluate(metrics, primary='thermal_robin', closest='thermal_robin_null')
    result = dict(observations=truth, metrics=metrics, refusals=refusals, decision=decision,
                  scope='Conditional nominal three-stake comparison, one glacier and year')
    with (receipt_path.parent / 'TRANSFER_RESULTS.json').open('x') as out:
        json.dump(result, out, indent=2)
    print(json.dumps(result, indent=2))
