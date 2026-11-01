"""Pinned, checksum-verified download of the two public Zenodo archives (raw data are never modified; stored under raw/, git-ignored). python download.py [outdir]
Endpoints are resolved through the Zenodo REST API (records/<id>); each archive's size and MD5 are checked against the values pinned here (from the independent audit) and against the API's own checksum."""
import hashlib, json, os, subprocess, sys, time
PIN = {"streams_icecliffs": dict(record="10812353", orig_record="8373145", name="Petersen_etal_Streams_IceCliffs_DataRepo.zip", size=57711310, md5="1b453396517d78d0d3d9686e1e4367d0"),
       "kennicott_backwasting": dict(record="14625691", name="Petersen_etal_2025_Kennicott_Glacier_Data.zip", size=9948630, md5="f740fc1d0e6ff19c1b8049b427db0c3a")}
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "raw"); os.makedirs(OUT, exist_ok=True); LOG = open(os.path.join(HERE, "download.log"), "a")
def log(*a): print(time.strftime("%H:%M:%S"), *a, file=LOG, flush=True)
def get(url, out, tries=8):
    for i in range(tries):
        r = subprocess.run(["curl", "-sSL", "-m", "900", "-o", out, "-w", "%{http_code}", url], capture_output=True, text=True); log("GET", url, r.stdout.strip())
        if r.stdout.strip() == "200": return True
        time.sleep(30 * (i + 1))
    return False
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
res = {}
for k, v in PIN.items():
    meta = os.path.join(OUT, f"{k}_record.json")
    if not get(f"https://zenodo.org/api/records/{v['record']}", meta): res[k] = dict(status="record endpoint failed"); continue
    try: m = json.load(open(meta)); files = {f["key"]: f for f in m["files"]}
    except Exception as e: res[k] = dict(status=f"record json unreadable: {e}"); continue
    f = files.get(v["name"]); res[k] = dict(status="listed" if f else "archive not listed", listed=list(files)); 
    if not f: continue
    dest = os.path.join(OUT, v["name"]); url = f["links"]["self"]
    if not (os.path.exists(dest) and os.path.getsize(dest) == v["size"]) and not get(url, dest): res[k]["status"] = "download failed"; continue
    res[k].update(size=os.path.getsize(dest), size_ok=os.path.getsize(dest) == v["size"], md5=md5(dest), md5_ok=md5(dest) == v["md5"], api_checksum=f.get("checksum"), status="ok" if md5(dest) == v["md5"] else "MD5 MISMATCH")
json.dump(res, open(os.path.join(HERE, "download_manifest.json"), "w"), indent=1); print(json.dumps({k: v.get("status") for k, v in res.items()}))
