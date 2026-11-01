"""Forcing/observation ingestion from the pinned archives (raw data never modified; no interpolation to fabricate measurements)."""
import csv, datetime as dt, io, zipfile, os
RAW = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_audit", "raw")
KENN = "Petersen_etal_2025_Kennicott_Glacier_Data.zip"; STREAMS = "Petersen_etal_Streams_IceCliffs_DataRepo.zip"
PINNED = {KENN: (9948630, "f740fc1d0e6ff19c1b8049b427db0c3a"), STREAMS: (57711310, "1b453396517d78d0d3d9686e1e4367d0")}
_VERIFIED = {}
def verify_pinned(name):
    """size and MD5 of a raw archive against the pinned values; cached per process (key: size, mtime). Raises on mismatch so a corrupted/substituted file is never ingested."""
    import hashlib
    p = os.path.join(RAW, name)
    if not os.path.exists(p): raise FileNotFoundError(f"{p} missing: run `python -m ice_cliff_survival data-download` first")
    st = os.stat(p); key = (p, st.st_size, st.st_mtime_ns)
    if key in _VERIFIED: return True
    size, md5 = PINNED[name]
    if st.st_size != size: raise ValueError(f"{name}: size {st.st_size} != pinned {size}")
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    if h.hexdigest() != md5: raise ValueError(f"{name}: MD5 {h.hexdigest()} != pinned {md5}; refusing to ingest")
    _VERIFIED[key] = True; return True
def open_zip(name):
    verify_pinned(name); return zipfile.ZipFile(os.path.join(RAW, name))
def load_backwasting(year=None, site=None):
    """list of dict(year, site, t_utc[list datetime], D_cm[list float], slope_deg, aspect_deg). D = emerged HORIZONTAL stake length [cm] (cumulative; ±2 cm, ±1 h). Blank rows dropped (reported by the audit)."""
    z = open_zip(KENN); rd = lambda n: z.read(n).decode("utf-8-sig", "replace"); N = [n for n in z.namelist() if "__MACOSX" not in n]
    sl = {r["Stake"]: r for r in csv.DictReader(io.StringIO(rd([n for n in N if n.endswith("ice_cliff_slope.csv")][0])))}; asp = {r["Stake"]: r for r in csv.DictReader(io.StringIO(rd([n for n in N if n.endswith("ice_cliff_aspect.csv")][0])))}; out = []
    for n in sorted(n for n in N if "Backwasting_Data" in n and n.endswith(".csv")):
        y = int(n.split("Backwasting_Data/")[1].split("/")[0]); s = n.rsplit("/", 1)[1][:-4]
        if (year is not None and y != year) or (site is not None and s != site): continue
        t, D = [], []
        for r in csv.DictReader(io.StringIO(rd(n))):
            if r["Datetime"].strip() and r["D"].strip(): t.append(dt.datetime.strptime(r["Datetime"].strip(), "%Y-%m-%d %H:%M:%S")); D.append(float(r["D"]))
        f = lambda tab: (float(tab[s][str(y)]) if tab.get(s, {}).get(str(y), "") not in ("", None) else None)
        out.append(dict(year=y, site=s, t_utc=t, D_cm=D, slope_deg=f(sl), aspect_deg=f(asp), sigma_D_cm=2.0, sigma_t_h=1.0))
    return out
def gauge_discharge(profiles_cm, length_m, times_s, surface_to_mean=1.0):
    """Q = (mean cross-section area) x (float velocity): profiles_cm = list of [(distance_cm, depth_cm), ...] per transect profile; times_s = float transit times over length_m. surface_to_mean converts surface-float speed to depth-mean (1.0 = none, as in the source)."""
    areas = []
    for p in profiles_cm:
        p = sorted(p); areas.append(sum(0.5 * (p[i][1] + p[i + 1][1]) * (p[i + 1][0] - p[i][0]) for i in range(len(p) - 1)) * 1e-4)
    A = sum(areas) / len(areas); v = sum(length_m / t for t in times_s) / len(times_s) * surface_to_mean; return dict(area_m2=A, velocity_m_s=v, Q_m3_s=A * v)
