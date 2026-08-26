from src.common.units import Q_
from src.stellar_params.stellar_data import StellarParams
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.chemical_env.chem_env_data import ChemEnvResult
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
    chemistry={
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
    local={
        "N2": 0.95,
        "CH4": 4.9e-2,
        "H2": 1.1e-3,
        "CO": 6.0e-5,
        "Ar": 4.32e-5,
        "C2H6": 1.1e-5,
        "C2H4": 5.0e-7,
        "C2H2": 3.0e-6,
        "C3H8": 6.0e-7,
        "HCN": 7.0e-7,
        "CH3CCH": 2.0e-8,
        "CHCCCH": 2.0e-8,
        "C6H6": 3.8e-9,
        "C2N2": 9.0e-10,
        "HCCCN": 4.4e-8,
        "H2O": 8.0e-9,
        "CO2": 1.3e-8
    },
    layered=None
)

def get_titan_stellar_params() -> StellarParams:
    return SOLAR_PARAMS

def get_titan_planetary_params() -> PlanetaryEnvironmentParams:
    return TITAN_PLANETARY_PARAMS

def get_titan_chem_env() -> ChemEnvResult:
    return TITAN_CHEM_ENV
