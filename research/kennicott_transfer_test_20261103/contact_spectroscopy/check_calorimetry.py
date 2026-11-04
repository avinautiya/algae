"""Check the temperature-derived Robin parameters against 2020 calorimetry."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'contact_model'))
import contact_experiment as model

def main(null=False):
    fit = json.loads((HERE / 'RESULTS_FEATURE_ONLY.json').read_text())['primary_DC_and_harmonics']['b0_null' if null else 'Robin_fit']
    plan = model.load_plan()
    site = plan['sites']['P_AWS_1']
    alpha = fit['alpha_m2_s']
    length = 0. if null else fit['b_kR_m']
    site['inverse_profile'] = 'two_sensor_2020b'
    plan['profiles']['two_sensor_2020b'] = {
        'alpha_layered_m2_s': [alpha, alpha], 'depths_m': [.05, .15, .22]
    }
    training = json.loads((ROOT / 'frozen_kennicott/TRAIN_RELEASED.json').read_text())
    target = training['rows'][0]['observed_net_ablation_mmwe']

    def evaluate(capacity):
        resistance = length / (alpha * capacity)
        result = model.forward(plan, 'P_AWS_1', model.PRIMARY, resistance, capacity)
        return dict(capacity_J_m3K=capacity, resistance_m2K_W=resistance,
                    melt_mmwe=result['physical_melt_mmwe'],
                    relative_energy_residual=result['relative_energy_residual'])

    low, high = .8e6, 3.5e6
    probes = [evaluate(c) for c in (low, 1.5e6, high)]
    calibrated = None
    if probes[0]['melt_mmwe'] <= target <= probes[-1]['melt_mmwe']:
        for _ in range(24):
            mid = evaluate((low + high) / 2)
            if abs(mid['melt_mmwe'] - target) <= .0036:
                calibrated = mid
                break
            if mid['melt_mmwe'] < target:
                low = mid['capacity_J_m3K']
            else:
                high = mid['capacity_J_m3K']
    return dict(scope='2020 development only, effective homogeneous Robin boundary',
                alpha_m2_s=alpha, extrapolation_length_m=length, target_mmwe=target,
                probes=probes, calibrated=calibrated, accessed_2023_melt=False)

if __name__ == '__main__':
    null = '--null' in sys.argv[1:]
    result = main(null)
    with (HERE / ('CALORIMETRY_NULL.json' if null else 'CALORIMETRY_CHECK.json')).open('x') as output:
        json.dump(result, output, indent=2)
    print(json.dumps(result, indent=2))
