"""SYNTHETIC demonstration of the reduced toe stock with ILLUSTRATIVE coefficients (not environmental, not fitted). Shows the regimes: export capacity above supply => stock stays below the covering volume; below => the covering volume is reached no later than the bound (Vcrit - Vs0)/(Smin - Cmax)."""
import json, sys
from .config import ToeConfig, Provenance
from .runner import ToeRun
from .closures import ConstantSources, ConstantCapacity
from .geometry import critical_solid_volume
from .scalar_model import accumulation_hitting_bound
def main(args):
    out = args[0] if args else "/tmp/ice_cliff_demo.json"; cfg = ToeConfig(W=10.0, porosity=0.3, H_cliff=15.0, theta_deg=60.0, beta_deg=35.0); Vc = critical_solid_volume(cfg.W, cfg.porosity, cfg.H_cliff, cfg.theta_deg, cfg.beta_deg); res = {}
    for name, S, C in (("export_dominated", 2e-5, 5e-5), ("supply_dominated", 6e-5, 2e-5)):
        run = ToeRun(cfg, ConstantSources({"S_mantle": S}), ConstantCapacity(C), dict(sources=Provenance("constant", "synthetic"), capacity=Provenance("constant", "synthetic")), 3600.0)
        run.advance(365 * 86400.0 * 2); res[name] = dict(Vcrit=Vc, final_Vs=run.stock.Vs, buried=run.history[-1]["fully_buried"], hit_time_bound_s=(None if accumulation_hitting_bound(0.0, Vc, S, C) == float("inf") else accumulation_hitting_bound(0.0, Vc, S, C)), first_buried_s=next((h["t"] for h in run.history if h["fully_buried"]), None), label="SYNTHETIC illustrative coefficients")
    json.dump(res, open(out, "w"), indent=1); print(json.dumps(res, indent=1)); return 0
