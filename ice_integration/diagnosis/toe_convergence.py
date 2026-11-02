"""Toe-contact convergence + wet/dry threshold sensitivity. Does NOT modify or overwrite existing demo runs.
  python -m ice_integration.diagnosis.toe_convergence settled OUT.json HDRY dx [tol]   -> settled-state f_toe on the INITIAL (identical) bed, P5 first-interval T
  python -m ice_integration.diagnosis.toe_convergence history OUT.json HDRY dx n       -> coupled history (reconstructed forcing) for n 6-h intervals, with per-interval settled flag
  python -m ice_integration.diagnosis.toe_convergence windows OUT.json                -> identical-window comparison of the EXISTING runs
HDRY is a module constant frozen into numba kernels (cache=True), so each value is run on a private copy of the reactivation package with its own name (cache cannot be stale)."""
import datetime as dt, importlib, json, os, re, shutil, sys, time
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.abspath(os.path.join(HERE,"..",".."))
def private_pkg(h):
    tag="rx_h%s"%re.sub(r"[^0-9a-z]","",("%g"%h)); dst=os.path.join(os.environ.get("TOE_SCRATCH","/tmp/toe_scratch"),tag)
    if not os.path.isdir(dst):
        os.makedirs(os.path.dirname(dst),exist_ok=True); shutil.copytree(os.path.join(ROOT,"ice_cliff_survival","reactivation"),dst,ignore=shutil.ignore_patterns("__pycache__","figures*","frozen"))
        p=os.path.join(dst,"swe2d.py"); s=open(p).read(); s=s.replace("HDRY = 1e-4","HDRY = %r"%h); open(p,"w").write(s)
    sys.path.insert(0,os.path.dirname(dst)); return importlib.import_module(tag+".cases"), importlib.import_module(tag+".swe2d")
def settled(out,h,dx,tol=2e-4):
    cases,sw=private_pkg(h); from ice_integration import forcing
    f=forcing.FieldForcing("P5"); a=f.integrate(dt.datetime(2019,7,20),dt.datetime(2019,7,20,6)); sim=cases.build("strip",dx=dx); sim.Tw=float(a.mean_temperature(True)); sim.Q=0.07
    t0=time.time(); r=sim.hydraulics(tmax=6000.0,tol=tol,cool=True); d=sim.diagnostics()
    res=dict(HDRY=h,dx=dx,tol=tol,Tw=sim.Tw,converged=bool(r["converged"]),hyd_time_s=r["hyd_time"],f_toe=d["f_toe"],Q_toe=d["Q_toe"],water_resid=d["water_resid"],energy_resid=d["energy_resid_rel"],wet_ice_area_toe_m2=d["wet_ice_area_toe"],wall_s=time.time()-t0)
    json.dump(res,open(out,"w"),indent=1,default=float); print(res)
def history(out,h,dx,n):
    cases,sw=private_pkg(h); from ice_integration import forcing, coupling
    f=forcing.FieldForcing("P5"); sim=cases.build("strip",dx=dx); drv=coupling.CoupledDriver(sim,dict(note="toe_convergence"),True,3000.0); sim.hydraulics(tmax=3000.0,cool=False); rows=[]; start=dt.datetime(2019,7,20)
    for k in range(n):
        a=f.integrate(start+dt.timedelta(hours=6*k),start+dt.timedelta(hours=6*(k+1))); t0=time.time(); d=drv.advance(a,0.07)
        rows.append(dict(k=k,f_toe=d["f_toe"],settled=sim.hyd_time<3000.0,hyd_time_s=sim.hyd_time,water_resid=d["water_resid"],wall_s=time.time()-t0)); json.dump(dict(HDRY=h,dx=dx,rows=rows,melt=drv.melt),open(out,"w"),indent=1,default=float)
    print(rows[-1])
def windows(out):
    R=os.path.join(ROOT,"ice_integration","results","demo"); res={}
    def cross(rows,th=1e-3):
        v=[r["f_toe"] for r in rows]
        for i in range(1,len(v)):
            if v[i]<th<=v[i-1]: return (i-1+(v[i-1]-th)/(v[i-1]-v[i]))*6.0
        return None
    for m in("strict","reconstructed"):
        a=json.load(open(os.path.join(R,"demo_dx0.05_%s.json"%m)))["rows"]; b=json.load(open(os.path.join(R,"demo_dx0.025_%s.json"%m)))["rows"]; res[m]={}
        for name,(i,j) in {"k0-5 (first 36 h)":(0,6),"k0-13 (active contact, 84 h)":(0,14),"k0-55 (all, diluted by zeros)":(0,56)}.items():
            ma=sum(r["f_toe"] for r in a[i:j])/(j-i); mb=sum(r["f_toe"] for r in b[i:j])/(j-i); res[m][name]=dict(mean_dx0p05=ma,mean_dx0p025=mb,abs_diff=ma-mb,rel_diff_to_fine=(ma-mb)/mb)
        res[m]["per_interval_k0-13"]=[dict(k=k,coarse=a[k]["f_toe"],fine=b[k]["f_toe"],abs=a[k]["f_toe"]-b[k]["f_toe"],rel=(a[k]["f_toe"]-b[k]["f_toe"])/b[k]["f_toe"] if b[k]["f_toe"]>1e-9 else None) for k in range(14)]
        ca,cb=cross(a),cross(b); res[m]["contact_loss_time_h_(f_toe<1e-3)"]=dict(dx0p05=ca,dx0p025=cb,abs_diff_h=None if None in(ca,cb) else ca-cb)
        res[m]["not_recorded"]="the demo rows did not store per-interval hydraulic convergence flags; settled-state evidence comes from the separate 'settled'/'history' runs"
    json.dump(res,open(out,"w"),indent=1,default=float)
    for m in res: print(m,{k:v for k,v in res[m].items() if k not in("per_interval_k0-13",)})
if __name__=="__main__":
    sys.path.insert(0,ROOT); c=sys.argv[1]
    if c=="settled": settled(sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[5]) if len(sys.argv)>5 else 2e-4)
    elif c=="history": history(sys.argv[2],float(sys.argv[3]),float(sys.argv[4]),int(sys.argv[5]))
    else: windows(sys.argv[2])
