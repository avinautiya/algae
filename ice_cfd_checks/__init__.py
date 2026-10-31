"""Independent analytical references and exported-data diagnostics; no CFD runtime."""

from .references import DamBreak, StefanReference, rectangular_jump

__all__ = ["DamBreak", "StefanReference", "rectangular_jump"]
