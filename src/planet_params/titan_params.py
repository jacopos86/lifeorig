from src.common.units import Q_
from src.stellar_params.stellar_data import StellarParams
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.chemical_env.chemical_environment import ChemEnvResult
from src.stellar_params.stellar_data import SOLAR_PARAMS

#
#   hardcoded Titan preset objects
#

TITAN_PLANETARY_PARAMS = PlanetaryEnvironmentParams(
    name="Titan",
    planet_radius=Q_(2574.73e3, "m"),
    planet_mass=Q_(1.3452e23, "kg"),
    orbital_distance=Q_(1.43353e12, "m"),
    rotation_period=Q_(15.945, "day"),
    chemical_env={
        "mode": "preset",
        "atomic_abundances": {},
        "chemical_species": []
    },
    obliquity=0.0049,
    eccentricity=0.0288,
    tidal_locked=True,
    day_night_contrast=0.0
)

TITAN_CHEM_ENV = ChemEnvResult(
    mode="preset",
    # TODO: ADD LOCAL CHEMICAL ENVIRONMENT FOR TITAN
    local=None,
    layered=None
)

def get_titan_stellar_params() -> StellarParams:
    return SOLAR_PARAMS

def get_titan_planetary_params() -> PlanetaryEnvironmentParams:
    return TITAN_PLANETARY_PARAMS

def get_titan_chem_env() -> ChemEnvResult:
    return TITAN_CHEM_ENV
