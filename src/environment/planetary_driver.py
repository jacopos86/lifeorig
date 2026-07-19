import numpy as np
from src.common.phys_constants import G
from src.chemical_env.chem_env_data import ChemEnvInput, ChemEnvResult
from src.chemical_env.chemical_environment import SimulateChemEnv
from src.utilities.logging_module import log
from src.planet_params.earth_params import get_earth_chem_env
from src.planet_params.titan_params import get_titan_chem_env
from src.planet_params.mars_params import get_mars_chem_env
from src.planet_params.europa_params import get_europa_chem_env
from src.planet_params.venus_params import get_venus_chem_env

#
#    PRESET chemistry
#

PRESET_CHEMISTRY = {
    "Titan": get_titan_chem_env,
    "Earth": get_earth_chem_env,
    "Mars": get_mars_chem_env,
    "Europa": get_europa_chem_env,
    "Venus": get_venus_chem_env
}

#
#   MAIN PLANETARY SOLVER DRIVER
#

def planetary_solver_driver(input_params):
    planet = input_params.planetary_data
    star = input_params.stellar_data
    # 1. Basic orbital/rotation derived quantities
    planet.surface_gravity = G * planet.planet_mass / planet.planet_radius**2
    planet.escape_velocity = (
        2.0 * G * planet.planet_mass / planet.planet_radius
    )**0.5
    planet.day_night_period = None if planet.tidal_locked else planet.rotation_period
    planet.seasonal_period = (
        2.0 * np.pi
        * (planet.orbital_distance**3 / (G * (star.mass + planet.planet_mass)))**0.5
    ).to("day")
    # 2. Build exochemistry input from planet object
    chem = planet.chemistry or {}
    mode = chem.get("mode")
    if mode == "preset":
        preset_loader = PRESET_CHEMISTRY.get(planet.name)
        if preset_loader is None:
            planet.chemical_stationary_config = ChemEnvResult(mode="preset")
        else:
            planet.chemical_stationary_config = preset_loader()
        return planet
    chem_input = ChemEnvInput(
        mode=mode,
        atomic_abundances=chem.get("atomic_abundances"),
        chemical_species=chem.get("chemical_species"),
        pressure=input_params.local_env_data.get("pressure"),
        temperature=input_params.local_env_data.get("temperature"),
    )
    if mode in {"local_equilibrium", "layered_equilibrium"}:
        planet.chemical_stationary_config = SimulateChemEnv(
            chem_input=chem_input,
            planet_data=planet,
            stellar_data=star,
            atmosphere_data=planet.atmosphere,
            output_dir=input_params.working_dir,
        ).run()
        return planet
    if mode == "layered_disequilibrium":
        log.error("layered_disequilibrium solver is not implemented")
    log.error(f"Unknown chemistry mode: {mode}")
