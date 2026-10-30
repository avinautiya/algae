"""Generate analytical requests or check solver exports without importing a solver."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from . import __version__
from .checks import (ENERGY_FIELDS, MASS_FIELDS, InvalidExport, audit_energy,
                     audit_mass, compare_rows)
from .references import analytical_cases


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def generate(destination):
    destination.mkdir(parents=True, exist_ok=True)
    cases = []
    for case in analytical_cases():
        path = destination / (case.name + ".reference.csv")
        rows = case.rows()
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        cases.append({"name": case.name, "description": case.description,
                      "parameters": case.parameters, "units": dict(time_s="s", **case.units),
                      "duration_s": case.duration_s, "samples": case.samples,
                      "reference_file": path.name,
                      "reference_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    write_json(destination / "manifest.json", {
        "version": __version__, "evidence_type": "analytical_verification_only",
        "not_observations": True, "cases": cases})
    print("Generated %d analytical cases in %s" % (len(cases), destination))


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise InvalidExport("Missing or duplicate CSV column names")
        rows = []
        for row in reader:
            if None in row or any(value is None or not value.strip() for value in row.values()):
                raise InvalidExport("Incomplete CSV row")
            rows.append({key: float(value) for key, value in row.items()})
        return rows


def evaluate(candidates, report, rtol, atol):
    if not candidates.is_dir():
        raise InvalidExport("Candidate directory does not exist")
    results = {}
    for case in analytical_cases():
        path = candidates / (case.name + ".csv")
        if not path.is_file():
            results[case.name] = {"passed": False, "error": "Missing solver export: " + path.name}
            continue
        try:
            actual = read_csv(path)
            result = compare_rows(case.rows(), actual, rtol=rtol, atol=atol)
            result["candidate_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            columns = set(actual[0])
            for label, fields, audit in (("mass_budget", MASS_FIELDS, audit_mass),
                                          ("energy_budget", ENERGY_FIELDS, audit_energy)):
                if columns.intersection(fields):
                    # melt_mass_kg is a required observable in heat_to_melt, not
                    # by itself an indication of a submitted global mass ledger.
                    requested = bool(columns.intersection(set(fields) - {"melt_mass_kg"}))
                    if requested:
                        result[label] = audit(actual)
                        result["passed"] = result["passed"] and result[label]["passed"]
                    else:
                        result[label] = {"status": "not_submitted"}
                else:
                    result[label] = {"status": "not_submitted"}
            results[case.name] = result
        except (ValueError, OSError) as error:
            results[case.name] = {"passed": False, "error": str(error)}
    summary = {"evidence_type": "analytical_verification_only",
               "passed": all(r["passed"] for r in results.values()), "cases": results}
    report.parent.mkdir(parents=True, exist_ok=True)
    write_json(report, summary)
    print("%d/%d cases passed; report: %s" %
          (sum(r["passed"] for r in results.values()), len(results), report))
    return 0 if summary["passed"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    gen = commands.add_parser("generate")
    gen.add_argument("destination", type=Path)
    check = commands.add_parser("evaluate")
    check.add_argument("candidates", type=Path)
    check.add_argument("--report", type=Path, required=True)
    check.add_argument("--rtol", type=float, default=0.01)
    check.add_argument("--atol", type=float, default=1e-8)
    args = parser.parse_args()
    if args.command == "generate":
        generate(args.destination)
        return 0
    try:
        return evaluate(args.candidates, args.report, args.rtol, args.atol)
    except (ValueError, OSError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
