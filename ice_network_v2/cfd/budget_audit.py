"""Exact discrete thermal budget of the scalarTransportFoam solution (steady, converged) for step_vs_ramp cases, reconstructed from the SAME quantities the solver used:
 equation per cell:  sum_f phi_f T_f  -  T_P sum_f phi_f  -  sum_f D snGrad_f |S_f| = 0     (div(phi,T) is a 'bounded' scheme: the T_P div(phi) term is part of the discretised equation)
 summing over cells, interior faces cancel pairwise; the global identity is
   R = sum_bnd [ phi_b T_b + D (T_P - T_b) delta_b |S_b| ]  -  sum_P T_P (sum_f phi_f)_P ,   phi_b positive OUT of the domain, delta_b = 1/max(n.d, 0.05|d|) (OpenFOAM nonOrthDeltaCoeffs on boundaries).
Terms reported (all per unit depth dz, nondimensional T_in = 1; denominators: P_in_adv = |phi_inlet| T_in):
 inlet advective, inlet CONDUCTIVE (axial conduction through the inlet face, nonzero at the heated-wall corner), outlet advective, outlet conductive (zero: zeroGradient), wall conduction (the melt-driving flux), lid (symmetry: 0), bounded term.
Also reports the earlier postprocessing ledger (Pin nominal, Pout from cell values u*T, wall flux with first-cell gradient) for comparison."""
import json, os, re, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import step_vs_ramp as S


def read_phi(path, patches, nfaces_int):
    t = open(path).read(); i = t.index("internalField"); j = t.index("(", i); head = t[i:j].split(); n = int(head[-1]); body = t[j + 1:].split("\n)")[0].strip().splitlines(); inner = np.array([float(x) for x in body[:n]])
    bf = t[t.index("boundaryField"):]; b = {}
    for name in patches:
        m = re.search(name + r"\s*\{(.*?)\n    \}", bf, re.S); blk = m.group(1)
        mm = re.search(r"value\s+nonuniform List<scalar>\s*(\d+)?\s*(?:\(([^)]*)\)|\n(\d+)\n\(([^)]*)\))", blk, re.S)
        mu = re.search(r"value\s+uniform\s+([-+0-9.eE]+);", blk)
        if mu and "nonuniform" not in blk: b[name] = float(mu.group(1)); continue
        mm = re.search(r"nonuniform List<scalar>\s*(\d+)\s*\(([^)]*)\)", blk, re.S) or re.search(r"nonuniform List<scalar>\s*\n(\d+)\n\(([^)]*)\)", blk, re.S)
        b[name] = np.array([float(x) for x in mm.group(2).split()]) if mm else None
    return inner, b


def audit(case, n, tl_flow="533"):
    pts, faces, owner, patches = S.mesh_data(case); nn, bb = S.read_list(f"{case}/constant/polyMesh/neighbour"); nbr = np.array([int(l) for l in bb.strip().splitlines()[:nn]])
    times = sorted([d for d in os.listdir(case) if d.replace(".", "").isdigit() and float(d) > 0], key=float); tT = times[-1]
    T = S.read_cell_field(f"{case}/{tT}/T", "sca"); ncell = len(T); C, fc = S.cell_centres(pts, faces, owner, nbr, ncell); dz = S.H / n; D = S.ALPHA
    inner, bphi = read_phi(f"{case}/{tT}/phi", ["inlet", "outlet", "walls", "lid"], nn)
    div = np.zeros(ncell); np.add.at(div, owner[:nn], inner); np.add.at(div, nbr, -inner)
    for name in ("inlet", "outlet"):
        s0, c0 = patches[name]; ph = bphi[name]
        for kk, k in enumerate(range(s0, s0 + c0)): div[owner[k]] += ph if np.isscalar(ph) else ph[kk]
    terms = dict(inlet_adv=0.0, inlet_cond=0.0, outlet_adv=0.0, outlet_cond=0.0, wall_cond=0.0, lid=0.0)
    def face_geom(k):
        p = pts[faces[k]]; z0 = [pp[:2] for pp in p if abs(pp[2]) < 1e-12]; e = z0[1] - z0[0]; L = np.linalg.norm(e); nrm = np.array([e[1], -e[0]]) / L; return L, nrm, fc[k][:2]
    for name, key in (("inlet", "inlet"), ("outlet", "outlet"), ("walls", "walls"), ("lid", "lid")):
        s, c = patches[name]; ph = bphi[name]
        for kk, k in enumerate(range(s, s + c)):
            L, nrm, Fc = face_geom(k); P = owner[k]; d = Fc - C[P][:2]
            if nrm @ d < 0: nrm = -nrm
            Sabs = L * dz; delta = 1.0 / max(nrm @ d, 0.05 * np.linalg.norm(d)); Tb = {"inlet": 1.0, "outlet": T[P], "walls": 0.0, "lid": T[P]}[name]
            phib = (ph if np.isscalar(ph) else ph[kk]) if ph is not None else 0.0
            if name == "inlet": terms["inlet_adv"] += phib * Tb; terms["inlet_cond"] += D * (T[P] - Tb) * delta * Sabs
            elif name == "outlet": terms["outlet_adv"] += phib * Tb; terms["outlet_cond"] += D * (T[P] - Tb) * delta * Sabs
            elif name == "walls": terms["wall_cond"] += D * (T[P] - Tb) * delta * Sabs
    bounded = float(np.sum(T * div)); R = sum(terms.values()) - bounded
    Pin = -terms["inlet_adv"]
    # previous postprocessing ledger
    r = S.analyse(case, n); Pin_old = r["Qin"] * 1.0; Pout_old = r["Qout"] * r["Tb_out"]
    # first-cell wall conduction used previously (same quantity but with n.d distance and cell value): recompute for comparison
    return dict(n=n, cells=ncell, time_T=tT, Pin_adv=Pin, terms=terms, bounded_term=bounded, R_abs=R, R_over_Pin=R / Pin, div_phi_max_rel=float(np.max(np.abs(div)) / (abs(bphi["inlet"]))),
                energy_into_wall_over_Pin=terms["wall_cond"] / Pin, inlet_cond_over_Pin=terms["inlet_cond"] / Pin, outlet_adv_over_Pin=terms["outlet_adv"] / Pin,
                closure_without_inlet_cond=(Pin - terms["outlet_adv"] - terms["wall_cond"]) / Pin, closure_with_inlet_cond=(Pin + terms["inlet_cond"] * -1 - terms["outlet_adv"] - terms["wall_cond"]) / Pin,
                old_ledger=(Pin_old - Pout_old - terms["wall_cond"]) / Pin_old)


if __name__ == "__main__":
    wd = sys.argv[1]; out = {}
    for kind in ("step", "ramp"):
        for n in (8, 16, 32):
            r = audit(f"{wd}/{kind}_n{n}", n); out[f"{kind}_n{n}"] = r
            print(kind, n, "R/Pin %.2e" % r["R_over_Pin"], "inletcond/Pin %.3e" % r["inlet_cond_over_Pin"], "wall/Pin %.5f" % r["energy_into_wall_over_Pin"], "outlet/Pin %.5f" % r["outlet_adv_over_Pin"], "closure(no inlet cond) %.2e" % r["closure_without_inlet_cond"], "bounded %.2e" % (r["bounded_term"] / r["Pin_adv"]), "divphi %.1e" % r["div_phi_max_rel"], flush=True)
    json.dump(out, open(os.path.join(HERE, "budget_audit.json"), "w"), indent=1, default=float)
