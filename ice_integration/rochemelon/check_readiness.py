"""Rochemelon independent-data readiness gate. Never scores; reports unresolved gates.
Usage: python -m ice_integration.rochemelon.check_readiness INPUT.json [--out report.json]
Exit 0 only if every gate is resolved (then the packet may be handed to predict_external)."""
import json, sys
from pathlib import Path
def _bad(v): return v is None or (isinstance(v,str) and (not v.strip() or v.startswith("REQUIRED")))
def _any_bad(x):
    if isinstance(x,dict): return any(_any_bad(v) for k,v in x.items() if not k.startswith("_"))
    if isinstance(x,list): return any(_any_bad(v) for v in x)
    return _bad(x)
def check(d):
    g=d.get("gates",{}); out={}
    def gate(n,ok,why): out[n]={"resolved":bool(ok),"reason":why}
    sb=g.get("section_binding",{}); gate("section_binding",not _any_bad(sb) and len(sb.get("predicted_section_names",[]))==3,"measured A/B/C must map to named prediction sections with distances")
    st=g.get("survey_timing",{}); gate("survey_timing",not _any_bad(st),"survey dates, hours and timing uncertainty required")
    it=g.get("inlet_thermal_path",{}); gate("inlet_thermal_path",not _any_bad(it),"inlet-to-target length and temperature source required")
    ge=g.get("geometry",{}); gate("geometry",not _any_bad(ge) and not ge.get("uses_rochemelon_fitted_coefficients",True),"independent geometry or qualified closure; no Rochemelon-fitted coefficients")
    al=d.get("alignment",{}); gate("alignment_flags",al.get("path_alignment_verified") is True and al.get("geometry_independent_of_bed_target") is True,"both alignment flags must be true")
    gate("records",bool(d.get("records")) and not _any_bad(d.get("records")),"forcing records with sections, no placeholders")
    ready=all(v["resolved"] for v in out.values())
    return {"ready":ready,"gates":out,"scored":False}
if __name__=="__main__":
    r=check(json.loads(Path(sys.argv[1]).read_text())); s=json.dumps(r,indent=1); print(s)
    if "--out" in sys.argv: Path(sys.argv[sys.argv.index("--out")+1]).write_text(s)
    sys.exit(0 if r["ready"] else 2)
