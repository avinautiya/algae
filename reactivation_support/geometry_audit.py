"""Measure realized pinned-snapshot geometry before interpreting refinement."""
from pathlib import Path
import json
from .snapshot import cases


def audit():
    rows=[]
    for dx in (.1,.05,.025,.0125):
        for layout in ('corridor','strip','porous_strip'):
            s=cases.build(layout,dx=dx)
            nominal_area=.5*s.Wl*s.Lr
            actual_area=s.deposit_cells*dx*dx
            rows.append(dict(layout=layout,dx=dx,deposit_cells=s.deposit_cells,
                nominal_covered_area_m2=nominal_area,actual_covered_area_m2=actual_area,
                cover_area_relative_error=actual_area/nominal_area-1,
                geometric_deposit_volume_m3=s.deposit_volume,
                grain_solid_volume_if_phi_applies_m3=s.deposit_volume*(1-.4) if layout=='porous_strip' else s.deposit_volume,
                computational_pool_width_m=s.ny*dx,
                represented_divider_width_m=dx,
                interpretation='porous geometric volume is bulk volume; solid grain volume requires porosity'))
    return rows

if __name__=='__main__':
    rows=audit();p=Path(__file__).parent/'records/geometry_audit.json';p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
