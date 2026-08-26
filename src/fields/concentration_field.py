from dataclasses import dataclass, field
import numpy as np
from src.common.units import Q_, ureg

@dataclass
class ConcentrationField:
    grid: object
    cell_volumes: Q_
    values: dict[str, Q_] = field(default_factory=dict)
    def __post_init__(self):
        self._validate_cell_volumes()
        for species, concentration in self.values.items():
            self._validate_species_concentration(species, concentration)
    # species
    def species(self):
        return set(self.values)
    # set species
    def set_species(self, species: str, concentration: Q_):
        self._validate_species_concentration(species, concentration)
        self.values[species] = concentration
    # get species
    def get_species(self, species: str, default=None):
        return self.values.get(species, default)
    # total number moles
    def total_moles(self, species: str):
        concentration = self.values.get(species)
        if concentration is None:
            return Q_(0.0, "mole")
        return np.sum(concentration * self.cell_volumes).to("mole")
    # total number molecules
    def total_molecules(self, species: str):
        return (self.total_moles(species) * ureg.avogadro_constant).to_base_units()
    # validation
    def _validate_cell_volumes(self):
        if not isinstance(self.cell_volumes, Q_):
            raise TypeError("cell_volumes must be a quantity")
        if self.cell_volumes.size != self.grid.n_cells:
            raise ValueError("cell_volumes must have one value per grid cell")
        if not self.cell_volumes.check("[length] ** 3"):
            raise ValueError("cell_volumes must have volume units")
    def _validate_species_concentration(self, species: str, concentration: Q_):
        if not isinstance(species, str):
            raise TypeError("species must be a string")
        if not isinstance(concentration, Q_):
            raise TypeError("concentration must be a quantity")
        if concentration.size != self.grid.n_cells:
            raise ValueError("concentration must have one value per grid cell")
        if not concentration.check("[substance] / [length] ** 3"):
            raise ValueError("concentration must have amount concentration units")