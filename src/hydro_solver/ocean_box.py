from dataclasses import dataclass, field


@dataclass
class OceanModel:
    """
    Conservative state of a well-mixed ocean.

    All inventories are defined per total planetary surface area.
    """
    # water mass per area (Kg/m^2)
    Sw: object
    # temperature
    temperature: object
    # fraction planet surface
    # covered by ocean f0 = A_ocean / A_planet
    f0: float = 1.0
    # dissolved chemicals -> amount chemical
    # per unit total planetary area
    dissolved_chem: dict[str, object] = field(
        default_factory=dict
    )
    def validate(self):
        if not 0.0 < self.f0 <= 1.0:
            log.error(
                "ocean_fraction must be in the interval (0, 1]"
            )
        water_inventory = self.Sw.to(
            "kilogram / meter**2"
        )
        if water_inventory.magnitude <= 0.0:
            log.error(
                "water_mass_per_area must be positive"
            )