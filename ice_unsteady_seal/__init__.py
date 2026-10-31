"""Finite-interval linear prototype of an unsteady ice-seal state model.

All variables are dimensionless. This is not a completed nonlinear breach model.
"""

from .model import BaseConfig, HydraulicMap, SealModel, build_base

__all__ = ["BaseConfig", "HydraulicMap", "SealModel", "build_base"]
