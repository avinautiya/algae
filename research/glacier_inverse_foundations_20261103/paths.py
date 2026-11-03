"""Locate the separately versioned Pirámide evidence package."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parent
DEFAULT=ROOT.parent/'piramide_thermal_test_20261103'
LOCAL=Path('/private/tmp/glacier-round2-20261103')
PIRAMIDE_ROOT=Path(os.environ.get('PIRAMIDE_RESEARCH_ROOT',str(DEFAULT if DEFAULT.is_dir() else LOCAL)))
def require_piramide():
    if not (PIRAMIDE_ROOT/'RELEASE_MANIFEST.json').is_file():
        raise FileNotFoundError('Set PIRAMIDE_RESEARCH_ROOT to the released Pirámide package or keep repository research folders as siblings')
    return PIRAMIDE_ROOT
