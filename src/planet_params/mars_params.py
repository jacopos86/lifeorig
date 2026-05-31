from src.common.units import Q_
from src.stellar_params.stellar_data import StellarParams
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.chemical_env.chemical_environment import ChemEnvResult
from src.stellar_params.stellar_data import SOLAR_PARAMS

#
#   hardcoded Mars preset objects
#

# Mars test mix, normalized to C = 1.0.
# This is a CO2-dominated thin atmosphere with N present as a minor component.

MARS_ATOMIC_ABUNDANCES = {
    "H": 6.3e-4,
    "He": 0.0,
    "C": 1.0,
    "N": 5.66e-2,
    "O": 2.002,
}

MARS_CHEMICAL_SPECIES = [
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

MARS_PLANETARY_PARAMS = PlanetaryEnvironmentParams(
    name="Mars",
    planet_radius=Q_(3.3895e6, "m"),
    planet_mass=Q_(6.4171e23, "kg"),
    orbital_distance=Q_(1.523679, "AU"),
    rotation_period=Q_(24.6229, "hour"),
    chemical_env={
        "mode": "layered_equilibrium",
        "atomic_abundances": MARS_ATOMIC_ABUNDANCES,
        "chemical_species": MARS_CHEMICAL_SPECIES,
    },
    obliquity=0.4396,
    eccentricity=0.0934,
    tidal_locked=False,
    day_night_contrast=0.0,
    atmosphere={
        "n_layers": 250,
        "z_max": Q_(120.0, "km"),
        "atmosphere_mass_fraction": 3.9e-8,
        "max_iter_loop": 100,
        "rel_tol": 1.0e-6,
        "abs_tol": 0.0,
        "damping_loop": 1.0,
    },
)

MARS_CHEM_ENV = ChemEnvResult(
    mode="layered_equilibrium",
    local=None,
    layered=None,
)

def get_mars_stellar_params() -> StellarParams:
    return SOLAR_PARAMS

def get_mars_planetary_params() -> PlanetaryEnvironmentParams:
    return MARS_PLANETARY_PARAMS

def get_mars_chem_env() -> ChemEnvResult:
    return MARS_CHEM_ENV
