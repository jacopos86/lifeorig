from dataclasses import dataclass
from src.common.units import Q_
from src.utilities.logging_module import log

#
#   data class for base planetary parameters
#

@dataclass
class PlanetaryEnvironmentParams:
    name: str
    planet_radius: Q_
    planet_mass: Q_
    orbital_distance: Q_
    rotation_period: Q_
    obliquity: float              # [rad]
    eccentricity: float           # dimensionless
    tidal_locked: bool = False
    day_night_contrast: float = 0.0   # 0 -> none, larger -> stronger contrast
    day_night_period: Q_ | None = None
    seasonal_period: Q_ | None = None
    # input configs
    atmosphere: dict | None = None
    hydro: dict | None = None
    chemistry: dict | None = None
    # planetary derived parameters
    surface_gravity: Q_ | None = None
    escape_velocity: Q_ | None = None
    # computed chemical equilibrium
    chemical_stationary_config: object | None = None
    # liquid water factor
    liquid_water_retention: float | None = None
    # log summary
    def log_summary(self):
        log.info("\t --- Planetary Parameters ---")
        log.info(f"\t planet_name        : {self.name}")
        log.info(f"\t planet_radius      : {self.planet_radius}")
        log.info(f"\t planet_mass        : {self.planet_mass}")
        log.info(f"\t orbital_distance   : {self.orbital_distance}")
        log.info(f"\t rotation_period    : {self.rotation_period}")
        log.info(f"\t obliquity          : {self.obliquity}")
        log.info(f"\t eccentricity       : {self.eccentricity}")
        log.info(f"\t tidal_locked       : {self.tidal_locked}")
