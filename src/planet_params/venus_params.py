from src.common.units import Q_
from src.stellar_params.stellar_data import StellarParams
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.chemical_env.chemical_environment import ChemEnvResult
from src.stellar_params.stellar_data import SOLAR_PARAMS

#
#   hardcoded Venus preset objects
#

# Venus test mix, normalized to C = 1.0.
# This is a dense CO2-dominated atmosphere with N present as a secondary
# component and trace H/S-free water chemistry available.

VENUS_ATOMIC_ABUNDANCES = {
    "H": 6.0e-5,
    "He": 0.0,
    "C": 1.0,
    "N": 7.2e-2,
    "O": 2.0,
}

VENUS_CHEMICAL_SPECIES = [
    "H",
    "H2",
    "He",
    "O",
    "O2",
    "OH",
    "H2O",
    "H2O2",
    "C",
    "CO",
    "CO2",
    "CH4",
    "N",
    "N2",
    "NH3",
    "NO",
    "NO2",
    "N2O",
    "HCN",
]

VENUS_PLANETARY_PARAMS = PlanetaryEnvironmentParams(
    name="Venus",
    planet_radius=Q_(6.0518e6, "m"),
    planet_mass=Q_(4.8675e24, "kg"),
    orbital_distance=Q_(0.723332, "AU"),
    rotation_period=Q_(243.025, "day"),
    chemical_env={
        "mode": "layered_equilibrium",
        "atomic_abundances": VENUS_ATOMIC_ABUNDANCES,
        "chemical_species": VENUS_CHEMICAL_SPECIES,
    },
    obliquity=3.096,
    eccentricity=0.0068,
    tidal_locked=False,
    day_night_contrast=0.0,
    atmosphere={
        "n_layers": 300,
        "z_max": Q_(250.0, "km"),
        "atmosphere_mass_fraction": 9.9e-5,
        "max_iter_loop": 100,
        "rel_tol": 1.0e-6,
        "abs_tol": 0.0,
        "damping_loop": 1.0,
    },
)

VENUS_CHEM_ENV = ChemEnvResult(
    mode="layered_equilibrium",
    local=None,
    layered=None,
)

def get_venus_stellar_params() -> StellarParams:
    return SOLAR_PARAMS

def get_venus_planetary_params() -> PlanetaryEnvironmentParams:
    return VENUS_PLANETARY_PARAMS

def get_venus_chem_env() -> ChemEnvResult:
    return VENUS_CHEM_ENV
