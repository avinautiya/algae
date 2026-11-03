"""JSON CLI: python -m thermal_memory INPUT.json OUTPUT.json"""
import argparse
import json
from pathlib import Path
from .solver import Column, SurfaceBoundary, ensemble


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    if args.output.exists():
        p.error('refusing to overwrite existing output')
    try:
        data = json.loads(args.input.read_text())
        boundary = SurfaceBoundary(**data['boundary'])
        columns = [Column(**x) for x in data['scenarios']]
        result = ensemble(columns, boundary, **data['options'])
        with args.output.open('x') as out:
            out.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    except (ValueError, TypeError, KeyError, OverflowError) as error:
        p.error(str(error))


if __name__ == '__main__':
    main()
