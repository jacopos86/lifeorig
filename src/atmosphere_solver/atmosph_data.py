import numpy as np
from dataclasses import dataclass, field
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

#
#   IR opacity settings
#

@dataclass
class IROpacitySettings:
    temperature_grid: object
    pressure_grid: object
    broadening: dict
    spectroscopy: dict
    cache_dir: str
    absorbers: tuple[str, ...]
    @classmethod
    def from_input(cls, data):
        return cls(
            temperature_grid=Q_(data["temperature_grid_K"], "K"),
            pressure_grid=Q_(data["pressure_grid_Pa"], "Pa"),
            broadening=dict(data["broadening"]),
            spectroscopy=dict(data["spectroscopy"]),
            cache_dir=data["cache_dir"],
            absorbers=tuple(data["absorbers"])
        )

#
#   Atmospheric surface (z=0) state
#

@dataclass(frozen=True)
class AtmosphereSurfaceState:
    temperature: object
    pressure: object
    mass_density: object
    wind_speed: object
    # chemical mole fractions
    mole_fractions: dict[str, object] = field(
        default_factory=dict
    )
    # radiation fluxes
    F_SW: object | None = None
    F_LW: object | None = None
    precipitation_flux: object | None = None
    # partial pressure
    def partial_pressure(self, species):
        try:
            mole_frac = self.mole_fractions[species]
        except KeyError:
            raise KeyError(
                f"No atmospheric mole fraction for {species}"
            ) from None
        return (
            mole_frac
            * self.pressure.to("pascal")
        ).to("pascal")