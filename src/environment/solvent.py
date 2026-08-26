from dataclasses import dataclass, field
from src.common.units import Q_

@dataclass
class SolventData:
    name: str
    composition: dict[str, float] = field(default_factory=dict)
    density: Q_ | None = None
    dynamic_viscosity: Q_ | None = None
    dielectric_constant: float | None = None
    diffusion_scale: float = 1.0
    polarity: float | None = None

#
#   Solvent class
#

class Solvent:
    def __init__(self, solvent_data: SolventData):
        self.name = solvent_data.name
        self.composition = dict(solvent_data.composition)
        self.density = solvent_data.density
        self.dynamic_viscosity = solvent_data.dynamic_viscosity
        self.dielectric_constant = solvent_data.dielectric_constant
        self.diffusion_scale = solvent_data.diffusion_scale
        self.polarity = solvent_data.polarity
    def species(self):
        return set(self.composition)
    # mole fractions
    def mole_fraction(self, species, default=0.0):
        return self.composition.get(species, default)
    # normalize composition
    def normalized_composition(self):
        total = sum(self.composition.values())
        if total <= 0.0:
            return {}
        return {
            species: fraction / total
            for species, fraction in self.composition.items()
        }
    def has_species(self, species):
        return species in self.composition
    def is_mixture(self):
        return len(self.composition) > 1
    def viscosity_scale(self):
        return self.dynamic_viscosity
    def diffusivity(self, base_diffusivity):
        return self.diffusion_scale * base_diffusivity