"""Nominal-value audit, not a replacement for authors' Monte Carlo products."""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def head_drop(z_up, h_up, z_down, h_down):
    if h_up < 0 or h_down < 0:
        raise ValueError('Depth must be nonnegative')
    return (z_up+h_up)-(z_down+h_down)

def f_nominal(drop, length, velocity, discharge, width=2.):
    if min(length, velocity, discharge, width) <= 0:
        raise ValueError('Positive hydraulic inputs required')
    area = discharge/velocity
    depth = area/width
    diameter = 4*area/(width+2*depth)
    return 2*9.81*drop*diameter/(length*velocity**2)

def f_stage_fixed(drop, length, velocity, discharge, h_up, h_down):
    depth = (h_up+h_down)/2
    if depth <= 0:
        raise ValueError('Positive mean depth required')
    return f_nominal(drop, length, velocity, discharge, width=(discharge/velocity)/depth)

def main():
    # Properties independent of the field observations.
    assert head_drop(5., .6, 4., .6) == 1.
    assert abs(head_drop(5., .6, 4., .4)-1.2) < 1e-12
    assert abs(head_drop(105., .6, 104., .4)-1.2) < 1e-12
    assert head_drop(4., .6, 4., .6) == 0.
    data = ROOT/'data/Channel_elevation'
    floors = {}
    for station in (3, 5):
        with (data/f'P{station}_channel_elevation.csv').open() as f:
            floors[station] = {r[0][:10]: float(r[1]) for r in list(csv.reader(f, delimiter=';'))[1:]}
    stages3 = {'2019-07-11': .23, '2019-07-16': .35, '2019-07-25': .55, '2019-07-30': .5, '2019-07-31': .38}
    stages5 = {'2019-07-11': .8, '2019-07-16': .6, '2019-07-25': .65, '2019-07-30': .6, '2019-07-31': .45}
    # Rounded observed velocities from source notebook cell87; exact travel times
    # and original stochastic particles are not reconstructed by this audit.
    velocity = [.387, .437, 1.057, 1.015, .583, .832, .564]
    q = [.206, .252, .811, 1.14, .923, .942, .619]
    with (ROOT/'products/Friction_factor_P5_P3.csv').open() as f:
        originals = list(csv.DictReader(f))
    rows = []
    for source, v, discharge in zip(originals, velocity, q):
        day = source['Date_Time'][:10]
        zu, zd = floors[5][day], floors[3][day]
        hu, hd = stages5[day], stages3[day]
        source_drop = zu+hu-zd+hd
        physical_drop = head_drop(zu, hu, zd, hd)
        assert abs(source_drop-physical_drop-2*hd) < 1e-10
        rows.append(dict(source_time=source['Date_Time'], source_exported_MC_mean_f=float(source['Friction_factor']),
                         z5=zu, z3=zd, h5=hu, h3=hd, Q_m3_s=discharge, velocity_rounded_m_s=v,
                         source_expression_drop_m=source_drop, piezometric_drop_m=physical_drop,
                         source_expression_nominal_f=f_nominal(source_drop, 253.7, v, discharge),
                         piezometric_nominal_f=f_nominal(physical_drop, 253.7, v, discharge),
                         source_expression_stage_fixed_nominal_f=f_stage_fixed(source_drop, 253.7, v, discharge, hu, hd),
                         piezometric_stage_fixed_nominal_f=f_stage_fixed(physical_drop, 253.7, v, discharge, hu, hd),
                         piezometric_over_source_head_ratio=physical_drop/source_drop,
                         nominal_roughness_ratio=math.sqrt(physical_drop/source_drop),
                         stage5_status='deduced, not measured' if day=='2019-07-11' else 'manual measured',
                         stage3_status='pressure sensor' if day=='2019-07-31' else 'manual measured'))
    result = dict(scope='Source algebra and nominal-value sensitivity only; not original Monte Carlo reproduction or validated corrected friction',
                  source_notebook_sha256='33750db4007b72418345f2ee3e93fb87165ac97eeb4a24e5d98fdd4dc149bf45',
                  source_cell_zero_based=131, property_checks=4, cases=rows)
    (ROOT/'HEAD_DROP_NOMINAL_AUDIT.json').write_text(json.dumps(result, indent=2)+'\n')
    with (ROOT/'HEAD_DROP_NOMINAL_CASES.csv').open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()
