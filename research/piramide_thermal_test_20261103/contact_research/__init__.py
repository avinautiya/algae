"""Additive contact-resistance foundation; no field validation or priority claim."""
from .contact_column import simulate
from .identifiability import invariants,steady_temperate_sensor_bound

__all__=['simulate','invariants','steady_temperate_sensor_bound']
