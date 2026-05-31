from dataclasses import dataclass
from src.common.units import Q_

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

#
#   chemical environment data results
#

@dataclass
class ChemEnvResult:
    mode: str
    local: object | None = None
    layered: object | None = None
    # get abundances
    def get_mole_fraction(self, species: str, default=None):
        if self.local == None:
            return default
        return self.local.get_mole_fraction(species, default)
    # get partial pressure
    def get_partial_pressure(self, species: str, default=None):
        if self.local == None:
            return default
        return self.local.get_partial_pressure(species, default)