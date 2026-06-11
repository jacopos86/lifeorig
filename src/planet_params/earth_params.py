from src.common.units import Q_
from src.stellar_params.stellar_data import StellarParams
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.chemical_env.chemical_environment import ChemEnvResult
from src.stellar_params.stellar_data import SOLAR_PARAMS

#
#   hardcoded Earth preset objects
#

# Primitive Earth test mix, normalized to N = 1.0.
# This intentionally avoids modern biogenic O2-rich air. It is a neutral to
# weakly reducing starting point with N2, CO2, H2O, H2, CO, CH4, and NH3
# available to the equilibrium chemistry.

EARTH_PRIMITIVE_ATOMIC_ABUNDANCES = {
    "H": 0.20,
    "He": 1.0e-8,
    "C": 0.03,
    "N": 1.0,
    "O": 0.078,
}

EARTH_CHEMICAL_SPECIES = [
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

# First-pass PHREEQC ocean composition preset. Values are dissolved elemental
# totals in mol/kgw, intended as a reasonable seawater-like starting point.
EARTH_OCEAN_COMPOSITION = {
    "Na": 0.47,
    "Cl": 0.55,
    "Mg": 0.053,
    "S": 0.028,
    "Ca": 0.010,
    "K": 0.010,
    "C": 0.002,
}

EARTH_PLANETARY_PARAMS = PlanetaryEnvironmentParams(
    name="Earth",
    planet_radius=Q_(6.371e6, "m"),
    planet_mass=Q_(5.972e24, "kg"),
    orbital_distance=Q_(1.0, "AU"),
    rotation_period=Q_(23.9344696, "hour"),
    chemical_env={
        "mode": "layered_equilibrium",
        "atomic_abundances": EARTH_PRIMITIVE_ATOMIC_ABUNDANCES,
        "chemical_species": EARTH_CHEMICAL_SPECIES,
    },
    obliquity=0.4091,
    eccentricity=0.0167,
    tidal_locked=False,
    day_night_contrast=0.0,
    atmosphere={
        "n_layers": 60,
        "z_max": Q_(100.0, "km"),
        "atmosphere_mass_fraction": 8.6e-7,
        "max_iter_loop": 100,
        # Convergence tolerances for pressure updates. Use a dimensionless
        # relative fraction and an absolute pressure tolerance Quantity.
        "rel_tol": 1.0e-4,
        "logp_tol": 1.0e-2,
        "abs_tol": Q_(1.0, "Pa"),
        "damping_loop": 0.3,
    },
    hydro={
        "ocean_composition": EARTH_OCEAN_COMPOSITION,
        "ocean_composition_units": "mol/kgw",
    },
)

EARTH_CHEM_ENV = ChemEnvResult(
    mode="layered_equilibrium",
    local=None,
    layered=None,
)

def get_earth_stellar_params() -> StellarParams:
    return SOLAR_PARAMS

def get_earth_planetary_params() -> PlanetaryEnvironmentParams:
    return EARTH_PLANETARY_PARAMS

def get_earth_chem_env() -> ChemEnvResult:
    return EARTH_CHEM_ENV
