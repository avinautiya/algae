"""Predeclared synthetic scale screen; these are not observed glacier inputs."""
import csv,json
from pathlib import Path
from storage_identity import Network,screening_bound
ROOT=Path(__file__).resolve().parent
cases=json.loads((ROOT/'campaign_cases.json').read_text());rows=[]
# Fixed explanatory footprint, not a field-derived geometry: 3m x0.5m face,
# lower5% coating fully wet, plus0.025m2 wet bare surface.
for r in cases:
    for h in [200.,2000.]:
        for excursion in [.5,10.]:
            area_debris=.075;area_bare=.025;area_total_face=1.5
            n=Network(r['capacity_J_m2_K']*area_debris,r['ice_conductance_W_m2_K']*area_debris,
                      h*area_bare,h*area_debris,1000.*4217.*r['inflow_m3_s'])
            row=dict(r,water_h_W_m2_K=h,assumed_terminal_excursion_K=excursion,
                     wetted_debris_area_m2=area_debris,wetted_bare_area_m2=area_bare,
                     total_face_footprint_m2=area_total_face,**n.diagnostics(),
                     **screening_bound(n,excursion,area_total_face))
            row['coating_footprint_equivalent_lowering_bound_m']=screening_bound(n,excursion,area_debris)['equivalent_lowering_bound_m']
            rows.append(row)
with (ROOT/'SCALE_SCREEN.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
result=dict(classification='assumed fixed-network scale bound only; not coupled outcomes or field data',
    cases=len(rows),maximum_face_equivalent_lowering_bound_m=max(r['equivalent_lowering_bound_m'] for r in rows),
    maximum_coating_equivalent_lowering_bound_m=max(r['coating_footprint_equivalent_lowering_bound_m'] for r in rows),
    largest_relaxation_time_s=max(r['relaxation_time_s'] for r in rows),
    caveat='Changing geometry invalidates the fixed-network identity; a larger irreversible effect needs demonstrated, numerically resolved amplification, not a claimed automatic bound on the nonlinear model.')
(ROOT/'SCALE_SCREEN.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
