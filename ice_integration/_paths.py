"""Explicit import bridge to the additive Codex packages (pinned 67f9b6ee…). Nothing is copied or modified."""
import os, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
for p in (os.path.join(ROOT, "ice_swiss_field_validation", "temperature_digitization"), os.path.join(ROOT, "ice_external_forecast", "model")):
    if p not in sys.path: sys.path.insert(0, p)
SWISS_SOURCE = os.path.join(ROOT, "ice_swiss_field_validation", "source")
