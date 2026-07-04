import numpy as np
from dataclasses import dataclass
from src.common.units import Q_

#
#   atmospheric physical data layer
#

@dataclass
class AtmLayerDyn:
    altitude: np.ndarray
    pressure: np.ndarray
    temperature: np.ndarray
    mean_molecular_mass: np.ndarray
    gravity: np.ndarray
    species_number_density: dict[str, np.ndarray]
    chemistry: object | None = None

#
#   atmosphere solver result
#

@dataclass
class AtmDynResult:
    atomic_abundances: dict[str, float]
    stellar_params: object
    uv_wavelength_grid: Q_
    ir_wavelength_grid: Q_
    stellar_B_lambda: Q_
    bond_albedo: float | None
    spectral_albedo: np.ndarray | None
    layers: AtmLayerDyn