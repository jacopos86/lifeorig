from src.common.units import Q_
from src.stellar_params.stellar_data import StellarParams
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.chemical_env.chemical_environment import ChemEnvResult
from src.stellar_params.stellar_data import SOLAR_PARAMS

#
#   hardcoded Europa preset objects
#

# Europa has a very tenuous O2-dominated exosphere sourced by surface ice
# radiolysis. This preset is mainly a low-atmosphere-mass stress test rather
# than a dense hydrostatic atmosphere benchmark.

EUROPA_ATOMIC_ABUNDANCES = {
    "H": 0.10,
    "He": 0.0,
    "C": 0.0,
    "N": 0.0,
    "O": 1.0,
}

EUROPA_CHEMICAL_SPECIES = [
    "H",
    "H2",
    "He",
    "O",
    "O2",
    "OH",
    "H2O",
    "H2O2",
]

# First-pass PHREEQC subsurface-ocean composition preset. Values are dissolved
# elemental totals in mol/kgw, intended as a simple sulfate/chloride icy-ocean
# starting point rather than a unique Europa composition.
EUROPA_OCEAN_COMPOSITION = {
    "Na": 0.05,
    "Cl": 0.05,
    "Mg": 0.01,
    "S": 0.02,
    "Ca": 0.002,
    "K": 0.001,
    "C": 0.001,
}

EUROPA_PLANETARY_PARAMS = PlanetaryEnvironmentParams(
    name="Europa",
    planet_radius=Q_(1.5608e6, "m"),
    planet_mass=Q_(4.7998e22, "kg"),
    orbital_distance=Q_(5.2044, "AU"),
    rotation_period=Q_(3.551, "day"),
    chemical_env={
        "mode": "layered_equilibrium",
        "atomic_abundances": EUROPA_ATOMIC_ABUNDANCES,
        "chemical_species": EUROPA_CHEMICAL_SPECIES,
    },
    obliquity=0.008,
    eccentricity=0.009,
    tidal_locked=True,
    day_night_contrast=0.0,
    atmosphere={
        "n_layers": 200,
        "z_max": Q_(300.0, "km"),
        "atmosphere_mass_fraction": 1.0e-16,
        "max_iter_loop": 100,
        "rel_tol": 1.0e-6,
        "abs_tol": 0.0,
        "damping_loop": 1.0,
    },
    hydro={
        "ocean_composition": EUROPA_OCEAN_COMPOSITION,
        "ocean_composition_units": "mol/kgw",
        "ocean_type": "subsurface",
    },
)

EUROPA_CHEM_ENV = ChemEnvResult(
    mode="layered_equilibrium",
    local=None,
    layered=None,
)

def get_europa_stellar_params() -> StellarParams:
    return SOLAR_PARAMS

def get_europa_planetary_params() -> PlanetaryEnvironmentParams:
    return EUROPA_PLANETARY_PARAMS

def get_europa_chem_env() -> ChemEnvResult:
    return EUROPA_CHEM_ENV
