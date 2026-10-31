"""A. Hydrostatic rest over a bump bed (stock interFoam, laminar, surface tension 0.07 N/m). 2-D, x in [0,1] m, bed y_b(x) = 0.05 exp(-((x-0.5)/0.1)^2), flat water level y = 0.2 m which is an exact MESH LINE (two blocks), atmosphere above.
Metrics at t = 1, 2 (and 5 s on the finest mesh): max |U| (spurious current), RMS of p_rgh/(rho g h_s) in the water, water-volume drift, mean free-surface height from the alpha column integral. Three resolutions."""
import json, os, shutil, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import common as C
L = 1.0; YS = 0.2; YT = 0.3
STEP_X = 0.5; STEP_H = 0.05          # orthogonal stepped bed: bed y = 0 for x < 0.5, y = 0.05 for x > 0.5 (Cartesian-aligned: no non-orthogonality). A smooth bump bed was tried first and FAILED (see FAILURES in the milestone report).
def make(case, ny_l, tend, maxco, wint=None):
    """ny_l = cells across the 0.2 m water layer (dy = 0.2/ny_l = dx); step height 0.05 and the atmosphere 0.1 are integer multiples of dy for ny_l multiple of 4."""
    shutil.rmtree(case, ignore_errors=True); dy = YS / ny_l; nxh = int(round(0.5 / dy)); dz = dy; n1 = int(round(STEP_H / dy)); n2 = ny_l - n1; n3 = int(round((YT - YS) / dy))
    V = [(0, 0), (.5, 0), (0, .05), (.5, .05), (1, .05), (0, .2), (.5, .2), (1, .2), (0, .3), (.5, .3), (1, .3)]; nv = len(V)
    verts = " ".join(f"({x} {y} 0)" for x, y in V) + " " + " ".join(f"({x} {y} {dz})" for x, y in V)
    blocks = [("A1", (0, 1, 3, 2), nxh, n1), ("A2", (2, 3, 6, 5), nxh, n2), ("B", (3, 4, 7, 6), nxh, n2), ("CA", (5, 6, 9, 8), nxh, n3), ("CB", (6, 7, 10, 9), nxh, n3)]
    def hexs(v): return f"hex ({v[0]} {v[1]} {v[2]} {v[3]} {v[0]+nv} {v[1]+nv} {v[2]+nv} {v[3]+nv})"
    def face(v, w):     # local faces of hex (l0..l3 bottom plane counter-clockwise, l4..l7 = +nv): w in bottom,top,left,right
        l = list(v) + [x + nv for x in v]
        return {"bottom": (l[0], l[1], l[5], l[4]), "top": (l[3], l[7], l[6], l[2]), "left": (l[0], l[4], l[7], l[3]), "right": (l[1], l[2], l[6], l[5])}[w]
    # hex local order (v0 v1 v2 v3) = (x0y0, x1y0, x1y1, x0y1): re-order the corner tuples above accordingly
    B = {n: (v[0], v[1], v[2], v[3]) for n, v, _, _ in blocks}
    f = lambda n, w: "(" + " ".join(map(str, face(B[n], w))) + ")"
    walls = [f("A1", "bottom"), f("A1", "left"), f("A1", "right"), f("A2", "left"), f("CA", "left"), f("B", "bottom"), f("B", "right"), f("CB", "right")]
    atm = [f("CA", "top"), f("CB", "top")]
    fb = " ".join("(" + " ".join(map(str, (v[0], v[3], v[2], v[1]))) + ") (" + " ".join(str(x + nv) for x in (v[0], v[1], v[2], v[3])) + ")" for _, v, _, _ in blocks)
    C.write(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"convertToMeters 1; vertices ({verts}); blocks ({' '.join(hexs(v) + f' ({nx} {ny} 1) simpleGrading (1 1 1)' for _, v, nx, ny in blocks)}); boundary (walls {{type wall; faces ({' '.join(walls)});}} atmosphere {{type patch; faces ({' '.join(atm)});}} frontAndBack {{type empty; faces ({fb});}});\n")
    C.physics(case); C.control(case, tend, 1e-4, maxco, wint or tend / 2)
    C.write(f"{case}/system/setFieldsDict", "dictionary", "setFieldsDict", f"defaultFieldValues (volScalarFieldValue alpha.water 0); regions (boxToCell {{box (-1 -1 -1) (2 {YS+1e-9} 1); fieldValues (volScalarFieldValue alpha.water 1);}});\n")
    C.write(f"{case}/0/alpha.water.orig", "volScalarField", "alpha.water", "dimensions [0 0 0 0 0 0 0]; internalField uniform 0; boundaryField { walls {type zeroGradient;} atmosphere {type inletOutlet; inletValue uniform 0; value uniform 0;} frontAndBack {type empty;} }\n")
    C.write(f"{case}/0/p_rgh", "volScalarField", "p_rgh", "dimensions [1 -1 -2 0 0 0 0]; internalField uniform 0; boundaryField { walls {type fixedFluxPressure; value uniform 0;} atmosphere {type totalPressure; p0 uniform 0; value uniform 0;} frontAndBack {type empty;} }\n")
    C.write(f"{case}/0/U", "volVectorField", "U", "dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0); boundaryField { walls {type noSlip;} atmosphere {type pressureInletOutletVelocity; value uniform (0 0 0);} frontAndBack {type empty;} }\n")
    return dict(dy=dy, nxh=nxh, n1=n1, n2=n2, n3=n3)
def read(path, kind):
    t = open(path).read(); i = t.index("internalField"); 
    if "uniform" in t[i:i + 40] and "nonuniform" not in t[i:i + 40]: return None
    j = t.index("(", i); n = int(t[i:j].split()[-1]); body = t[j + 1:].split("\n)")[0].strip().splitlines()
    return np.array([[float(v) for v in l.strip("()").split()] for l in body[:n]]) if kind == "vec" else np.array([float(l) for l in body[:n]])
def analyse(case, g, times):
    """cell order: blocks A1, A2, B, CA, CB each x-fastest. volumes uniform dy^3. Water cells: A1, A2, B (y < 0.2)."""
    dy = g["dy"]; nxh = g["nxh"]; V = dy ** 3; nA1 = nxh * g["n1"]; nA2 = nxh * g["n2"]; nB = nA2; nw = nA1 + nA2 + nB; out = {}
    vol0 = nw * V
    for t in times:
        a = read(f"{case}/{t}/alpha.water", "sca"); u = read(f"{case}/{t}/U", "vec"); p = read(f"{case}/{t}/p_rgh", "sca")
        assert np.all(np.isfinite(a)) and np.all(np.isfinite(u)) and np.all(np.isfinite(p)), "NaN in fields"
        umag = np.linalg.norm(u[:, :2], axis=1); rgh = 1000 * 9.81 * YS
        out[str(t)] = dict(max_U=float(umag.max()), max_U_water=float(umag[:nw].max()), water_volume_rel_drift=float((a.sum() * V - vol0) / vol0), alpha_min=float(a.min()), alpha_max=float(a.max()),
                           p_rgh_std_over_rho_g_hs=float(p[:nw].std() / rgh), p_rgh_mean_over_rho_g_hs=float(p[:nw].mean() / rgh), p_rgh_dev_max_over_rho_g_hs=float(np.max(np.abs(p[:nw] - rgh)) / rgh))
    return out
if __name__ == "__main__":
    wd = sys.argv[1]; res = {}
    for nyl in (16, 32, 64):
        case = f"{wd}/hs_{nyl}"; g = make(case, nyl, 2.0, 0.2, 1.0); C.run(["blockMesh"], case, "bm.log")
        shutil.copy(f"{case}/0/alpha.water.orig", f"{case}/0/alpha.water"); C.run(["setFields"], case, "sf.log"); C.run(["checkMesh"], case, "cm.log"); C.run(["interFoam"], case, "run.log")
        res[f"ny{nyl}"] = analyse(case, g, [1, 2]); res[f"ny{nyl}"]["dx"] = g["dy"]; print(nyl, json.dumps(res[f"ny{nyl}"]), flush=True)
        json.dump(res, open(os.path.join(HERE, "hydrostatic_results.json"), "w"), indent=1)
