"""Read immutable PROMBIO 2024 records; no invented coordinates or pseudocounts."""
from __future__ import annotations
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path


def observation(value):
    """Preserve missing and detection-limit observations distinctly from exact zero."""
    raw = value.strip()
    if raw.lower() in {"", "na", "n/a", "nan"}:
        return dict(kind="missing", value=None, upper_bound=None, raw=value)
    if raw.lower() == "no counts":
        return dict(kind="unresolved_annotation", value=None, upper_bound=None, raw=value,
                    note="do not infer zero versus uncounted without the protocol")
    if raw.startswith("<"):
        upper = float(raw[1:].replace(",", "."))
        if not math.isfinite(upper) or upper <= 0:
            raise ValueError("invalid detection limit")
        return dict(kind="left_censored", value=None, upper_bound=upper, raw=value)
    number = float(raw.replace(",", "."))
    if not math.isfinite(number) or number < 0:
        raise ValueError("finite nonnegative counts required")
    return dict(kind="observed", value=number, upper_bound=None, raw=value)


def load_2024(directory=None):
    directory = Path(directory) if directory else Path(__file__).parent / "data"
    path = directory / "prombio_2024_original.csv"
    metadata = json.loads((directory / "prombio_2024_metadata.json").read_text())
    version = metadata["data"]["latestVersion"]
    source = next(f["dataFile"] for f in version["files"] if f["dataFile"]["filename"].endswith(".csv"))
    if hashlib.md5(path.read_bytes()).hexdigest() != source["md5"]:
        raise ValueError("provider checksum mismatch")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter=";"))
    headers = rows[1]
    records = []
    for line, row in enumerate(rows[2:], start=3):
        if not any(row):
            continue
        if len(row) != len(headers):
            raise ValueError(f"source row {line} has unexpected column count")
        raw = dict(zip(headers, row))
        record = dict(sample_id=raw["Sample unique ID"], source_line=line,
                      station_area=raw["Sample"].strip().replace("-", "_"),
                      date=datetime.strptime(raw["date"], "%d-%m-%Y").date().isoformat(),
                      sampling_stratum=raw["Randomized grid or under radiometer"].strip(),
                      ice_cells_ml=observation(raw["ice algae (cells/ml)"]),
                      snow_cells_ml=observation(raw["snow algae (cells/ml)"]),
                      total_cells_ml=observation(raw["ice+snow algae (cells/ml)"]),
                      latitude=None, longitude=None, position_status="not supplied in this table",
                      spectral_observation_status="not supplied in this table",
                      raw=raw)
        records.append(record)
    if len({r["sample_id"] for r in records}) != len(records):
        raise ValueError("duplicate sample IDs")
    return records, dict(doi=metadata["data"]["persistentUrl"], version=f"{version['versionNumber']}.{version['versionMinorNumber']}",
                         file_id=source["id"], sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         license=version["license"]["name"], units="cells per mL; verify sampling/meltwater convention before model conversion")


if __name__ == "__main__":
    from collections import Counter
    records, provenance = load_2024()
    print(json.dumps(dict(provenance=provenance, n_records=len(records),
                          stations=dict(Counter(r["station_area"] for r in records)),
                          strata=dict(Counter(r["sampling_stratum"] for r in records))), indent=2))
