"""Form three candidate predictions from explicitly supplied temperature-only fits."""
import argparse
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'contact_model'))
import contact_experiment as model

def predict(fits, capacity):
    sites = ('P11', 'P12', 'P13')
    if set(fits) != set(sites) or not math.isfinite(capacity) or not .8e6 <= capacity <= 3.5e6:
        raise ValueError('Require all three source-bound fits and an admissible training capacity')
    plan = model.load_plan()
    predictions = []
    for site in sites:
        fit = fits[site]
        alpha = fit['alpha_m2_s']
        length = fit['b_kR_m']
        if not math.isfinite(alpha) or not 1e-8 <= alpha <= 5e-6 or not math.isfinite(length) or not 0 <= length <= .3:
            raise ValueError(f'{site}: parameters outside the declared fit bounds')
        resistance = length / (alpha * capacity)
        profile = f'robin_fit_{site}'
        original_profile = plan['sites'][site]['inverse_profile']
        plan['sites'][site]['inverse_profile'] = profile
        plan['profiles'][profile] = {
            'alpha_layered_m2_s': [alpha, alpha],
            'depths_m': plan['profiles'][original_profile]['depths_m'],
        }
        result = model.forward(plan, site, model.PRIMARY, resistance, capacity)
        predictions.append(dict(site=site, start='2023-06-08', end='2023-09-02',
                                predicted_mmwe=result['physical_melt_mmwe'],
                                alpha_m2_s=alpha, b_kR_m=length, resistance_m2K_W=resistance,
                                relative_energy_residual=result['relative_energy_residual']))
    return predictions

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fits', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--calibration', type=Path, default=HERE / 'CALORIMETRY_CHECK.json')
    args = parser.parse_args()
    calibration = json.loads(args.calibration.read_text())['calibrated']
    if calibration is None:
        raise ValueError('2020 capacity calibration refused')
    rows = predict(json.loads(args.fits.read_text()), calibration['capacity_J_m3K'])
    with args.output.open('x') as output:
        json.dump(dict(scope='Candidate predictions before 2023 outcome inspection', rows=rows), output, indent=2)
