"""Sibling-package resolution without numerical/source mutations."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def round3_root():
    candidates=[Path(os.environ['GLACIER_ROUND3_ROOT'])]if os.environ.get('GLACIER_ROUND3_ROOT')else []
    candidates += [ROOT.parent/'glacier_inverse_foundations_20261103',ROOT.parent/'glacier-round3-20261103']
    for p in candidates:
        if (p/'RELEASE_MANIFEST.json').is_file():return p
    raise FileNotFoundError('Keep the glacier_inverse_foundations_20261103 sibling directory or set GLACIER_ROUND3_ROOT')

def round2_root():
    candidates=[Path(os.environ['PIRAMIDE_RESEARCH_ROOT'])]if os.environ.get('PIRAMIDE_RESEARCH_ROOT')else []
    candidates += [ROOT.parent/'piramide_thermal_test_20261103',ROOT.parent/'glacier-round2-20261103']
    for p in candidates:
        if (p/'RELEASE_MANIFEST.json').is_file():return p
    raise FileNotFoundError('Keep the piramide_thermal_test_20261103 sibling directory or set PIRAMIDE_RESEARCH_ROOT')
