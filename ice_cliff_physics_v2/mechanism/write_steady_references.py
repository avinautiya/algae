"""Write analytic forecasts before independent resolved-slab comparisons."""
import csv,hashlib,json
from pathlib import Path
from wet_resistance import prediction
ROOT=Path(__file__).resolve().parent
freeze=json.loads((ROOT/'WET_BOUNDARY_FREEZE.json').read_text())
for n,key in [('WET_RESISTANCE_SPEC.md','spec_sha256'),('wet_boundary_cases.json','cases_sha256')]:
    if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=freeze[key]:raise ValueError('frozen input changed')
rows=[]
for c in json.loads((ROOT/'wet_boundary_cases.json').read_text()):
    v=prediction(water_capacity_W_K=c['water_capacity_W_K'],film_conductance_W_K=c['film_conductance_W_K'],ice_conductance_W_K=c['ice_conductance_W_K'],inlet_superheat_K=c['inlet_superheat_K'],absorbed_external_heat_W=c['absorbed_source_W'])
    A=c['total_covered_area_m2'];Aw=c['wet_area_m2'];qdry=c['absorbed_source_W']/Aw*(A-Aw)
    rows.append(dict(c,**v,wet_footprint_ice_flux_W_m2=v['covered_ice_power_W']/Aw,
                     whole_covered_ice_flux_W_m2=(v['covered_ice_power_W']+qdry)/A,
                     whole_covered_bare_proxy_flux_W_m2=(v['bare_ice_power_W']+qdry)/A,
                     whole_covered_dry_port_flux_W_m2=c['absorbed_source_W']/Aw,
                     outgoing_water_warming_K=v['outgoing_water_temperature_C']-c['inlet_superheat_K']))
with (ROOT/'ANALYTIC_STEADY_REFERENCES.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(ROOT/'ANALYTIC_REFERENCE_NOTICE.json').write_text(json.dumps(dict(classification='analytic constitutive forecasts, not observations or resolved-solver outcomes',cases=len(rows),maximum_energy_residual_W=max(abs(r['energy_residual_W']) for r in rows),expected_warming_cases=sum(r['wet_minus_dry_ice_power_W']>0 for r in rows),expected_cooling_cases=sum(r['wet_minus_dry_ice_power_W']<0 for r in rows)),indent=2)+'\n')
print('Wrote',len(rows),'analytic reference cases')
