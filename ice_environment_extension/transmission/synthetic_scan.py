"""Synthetic sensitivity only: no calibrated aperture range or discovery claim."""
import json
import math
from pathlib import Path
from transmission import Route, solve_common_stage


def run_scan():
    rho_cp = 1000.0 * 4217.0  # J/(m³ K), assumed constant liquid water
    ua, tw, inflow = 200.0, 0.5, 0.07
    rows = []
    for aperture in [0.0] + [10**(-5 + i/10) for i in range(41)]:
        routes = [Route('toe', 0.05, 0.2, 0.5, 1.7, aperture, 0.6),
                  Route('bypass', 0.0, 0.0, 0.5, 1.7)]
        solution = solve_common_stage(routes, inflow, 0, flow_tolerance_m3_s=1e-12)
        q = solution.route_flows_m3_s['toe']
        capacity = rho_cp*q
        inlet_heat = capacity*tw
        # Steady plug flow past an isothermal toe: T_out/T_in=exp(-UA/(rho cp Q)).
        residual_heat = inlet_heat*math.exp(-ua/capacity) if capacity else 0.0
        wetted=min(aperture,max(solution.stage_m-.05,0.))
        speed=q/(.5*wetted) if wetted>0 else None
        diameter=2*.5*wetted/(.5+wetted) if wetted>0 else None
        rows.append(dict(aperture_m=aperture, stage_m=solution.stage_m,
                         toe_flow_m3_s=q, bypass_flow_m3_s=solution.route_flows_m3_s['bypass'],
                         mass_residual_m3_s=solution.residual_m3_s,
                         toe_inlet_advected_heat_W=inlet_heat,
                         toe_floor_heat_W=inlet_heat*(-math.expm1(-ua/capacity)) if capacity else 0.0,
                         remaining_advected_heat_W=residual_heat,
                         opening_mean_speed_m_s=speed,
                         opening_Re_Dh_indicator=speed*diameter/1.79e-6 if speed is not None else None))
    return dict(label='SYNTHETIC sensitivity scan; NOT a defensible site aperture range',
                closure='Assumed serial plug-flow toe-floor UA; toe_floor_heat_W is actual modeled floor uptake, remaining heat is available downstream',
                applicability_warning='Thin-opening hydrostatic rating with constant Cd. Finite under-rock path length, viscous friction, roughness and turbulence are not resolved. Especially tiny-aperture results are NOT physical under-rock predictions. Reynolds indicator is diagnostic only.',
                constants=dict(inflow_m3_s=inflow, toe_UA_W_K=ua, inlet_excess_temperature_K=tw,
                               rho_cp_J_m3_K=rho_cp, downstream_stage_m=0,
                               route_widths_m=[0.5, 0.5], toe_floor_m=0.05,
                               bypass_floor_m=0, toe_debris_height_m=0.2,
                               Cw_sqrt_m_per_s=1.7, Cd=0.6), rows=rows)


def write_svg(data, path):
    rows = data['rows'][1:]
    ymax = max(r['toe_floor_heat_W'] for r in rows)
    points = ' '.join(f"{80+600*(math.log10(r['aperture_m'])+5)/4:.2f},{330-240*r['toe_floor_heat_W']/ymax:.2f}" for r in rows)
    ticks = ''.join(f'<text x="{80+150*i}" y="352" text-anchor="middle">10^{i-5}</text>' for i in range(5))
    path.write_text(f'''<svg xmlns="http://www.w3.org/2000/svg" width="780" height="450" viewBox="0 0 780 450">
<rect width="780" height="450" fill="white"/>
<g font-family="sans-serif" fill="#172b3a"><text x="30" y="28" font-size="19">SYNTHETIC aperture sensitivity</text>
<text x="30" y="51" font-size="13">Assumed serial UA=200 W/K; NOT a defensible site range</text>
<path d="M80 80 V330 H680" fill="none" stroke="#172b3a"/>
<polyline points="{points}" fill="none" stroke="#007e87" stroke-width="3"/>
<text x="70" y="334" text-anchor="end">0</text><text x="70" y="94" text-anchor="end">{ymax:.0f}</text>
{ticks}<text x="380" y="384" text-anchor="middle">Toe aperture height (m), logarithmic axis</text>
<text transform="translate(22 220) rotate(-90)" text-anchor="middle">Modeled toe-floor heat uptake (W)</text>
<text x="80" y="414" font-size="13">Zero aperture: Q_toe=0, heat=0 (omitted from logarithmic axis).</text>
<text x="80" y="434" font-size="12">Thin-opening model only; finite-path viscous resistance omitted. Tiny-gap values unvalidated.</text></g></svg>''')


if __name__ == '__main__':
    data = run_scan()
    root = Path(__file__).resolve().parent
    (root/'synthetic_scan.json').write_text(json.dumps(data, indent=2)+'\n')
    write_svg(data, root/'synthetic_scan.svg')
    print(json.dumps({'samples': len(data['rows']), 'zero': data['rows'][0],
                      'smallest_positive': data['rows'][1], 'largest': data['rows'][-1]}, indent=2))
