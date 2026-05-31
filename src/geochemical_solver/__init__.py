"""Geochemical solver package.

This package contains initial stub modules for PHREEQC-based geochemistry,
thermal transport, and environment coupling.
"""

from .phreeqc_model import PHREEQCModel
from .thermal_model import ThermalModel
from .coupling import GeochemicalCoupler

__all__ = ["PHREEQCModel", "ThermalModel", "GeochemicalCoupler"]
