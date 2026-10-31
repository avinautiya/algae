"""Rate-limit-tolerant retrieval of the Ogier et al. ETH Research Collection record (DSpace REST). Logs every endpoint and HTTP status; downloads README first, then the archive. Writes into ./raw (untrusted data: never executed)."""
import json, os, subprocess, sys, time
H = os.path.dirname(os.path.abspath(__file__)); RAW = os.path.join(H, "raw"); os.makedirs(RAW, exist_ok=True); LOG = open(os.path.join(H, "fetch.log"), "a")
def log(*a): print(time.strftime("%H:%M:%S"), *a, file=LOG, flush=True)
def get(url, out, tries=12, wait=60):
    for i in range(tries):
        p = subprocess.run(["curl", "-sS", "-L", "-m", "900", "-o", out, "-w", "%{http_code}", url], capture_output=True, text=True); code = p.stdout.strip()
        log("GET", url, "->", code, os.path.getsize(out) if os.path.exists(out) else -1)
        if code == "200": return True
        time.sleep(wait * (1 + i // 3))
    return False
base = "https://www.research-collection.ethz.ch/server/api/core/items/27f98ac0-2289-4539-9fb6-c6ee6308600a"
if get(base + "/bundles", f"{RAW}/bundles.json"):
    for b in json.load(open(f"{RAW}/bundles.json"))["_embedded"]["bundles"]:
        if b["name"] != "ORIGINAL": continue
        if get(b["_links"]["bitstreams"]["href"], f"{RAW}/bitstreams.json"):
            for bs in json.load(open(f"{RAW}/bitstreams.json"))["_embedded"]["bitstreams"]:
                log("bitstream", bs["name"], bs.get("sizeBytes"), bs["uuid"], bs.get("checkSum"))
                get(bs["_links"]["content"]["href"], f"{RAW}/{bs['name']}", tries=6, wait=120)
