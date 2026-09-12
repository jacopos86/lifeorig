from dataclasses import dataclass
import numpy as np
from src.common.units import Q_
from src.utilities.logging_module import log

#
#    data input for environment calculations
#

@dataclass
class ChemEnvInput:
    mode: str
    atomic_abundances: dict[str, float]
    chemical_species: list[str] | None = None
    pressure: Q_ | None = None
    temperature: Q_ | None = None
    reaction_source_files: list | None = None

#
#   chemical environment data results
#

@dataclass
class ChemEnvResult:
    mode: str
    local: object | None = None
    layered: object | None = None
    # get composition at a location
    def atmosph_composition(self, z=None, default=None):
        if self.local is not None:
            return self._local_composition(default=default)
        if self.layered is not None:
            return self._layered_composition(z=z, default=default)
        return default
    def composition_at(self, reservoir: str = "atmosphere", z=None, default=None):
        if reservoir != "atmosphere":
            return default
        return self.atmosph_composition(z=z, default=default)
    # LOCAL COMPOSITION 
    def _local_composition(self, default=None):
        if isinstance(self.local, dict):
            return dict(self.local)
        return dict(self.local.mole_fractions)
    # LAYERED COMPOSITION
    def _layered_composition(self, z=None, default=None):
        layers = self.layered.layers
        profiles = layers.chemistry.mole_fraction_profiles
        altitude = layers.altitude
        if not isinstance(altitude, Q_):
            log.error("layer altitude must be a quantity")
        if altitude.size == 0:
            return default
        if z is None:
            layer_index = int(np.argmin(altitude.to_base_units().magnitude))
        else:
            if not isinstance(z, Q_):
                log.error("z must be None or a quantity with altitude units")
            layer_distance = abs(altitude - z.to(altitude.units))
            layer_index = int(np.argmin(layer_distance.to_base_units().magnitude))
        return {
            species: float(np.asarray(profile, dtype=float)[layer_index])
            for species, profile in profiles.items()
        }
