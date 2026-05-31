from src.utilities.logging_module import log
from src.common.phys_constants import G
from src.planet_params.planetary_params import DerivedPlanetQuantities, PlanetaryEnvironmentParams
from src.chemical_env.chemical_environment import ChemEnvResult
from src.environment.liquid_level_model import derive_liquid_base_factor, derive_liquid_amplitude_factor

#
#   derive planetary parameters
#

def derive_planet_env_data(
        env_model: str,
        env_input_params: dict,
        planet_data: PlanetaryEnvironmentParams,
        chem_env_params: ChemEnvResult
    ) -> DerivedPlanetQuantities:
    # compute gravity acc.
    g = G * planet_data.planet_mass / planet_data.planet_radius**2
    # escape velocity
    vesc = (2.0 * G * planet_data.planet_mass / planet_data.planet_radius)**0.5
    # water base factor evaluation
    liquid_level_params = env_input_params.get("solvent_data").liquid_level_params
    water_base_factor = liquid_level_params.base_level
    if water_base_factor is None:
        water_base_factor = derive_liquid_base_factor(
            env_model=env_model,
            env_data=env_input_params,
            planet_data=planet_data,
            chemical_env=chem_env_params,
            vesc=vesc
        )
        log.info(f"\t water base factor: {water_base_factor}")
    # water amplitude factor
    water_amplitude_factor = liquid_level_params.amplitude
    if water_amplitude_factor is None:
        water_amplitude_factor = derive_liquid_amplitude_factor(
            env_model=env_model,
            planet_data=planet_data,
            water_base_factor=water_base_factor
        )
        log.info(f"\t water amplitude factor: {water_amplitude_factor}")
    radiation_amplitude_factor = (
        1.0 + 3.0 * planet_data.eccentricity
    ) * planet_data.day_night_contrast
    # atmosphere retention proxy for shielding
    radiation_base_factor = 1.0 / water_base_factor if water_base_factor else 1.0

    return DerivedPlanetQuantities(
        gravity=g,
        escape_velocity=vesc,
        water_base_factor=water_base_factor,
        water_amplitude_factor=water_amplitude_factor,
        radiation_base_factor=radiation_base_factor,
        radiation_amplitude_factor=radiation_amplitude_factor,
        day_night_period=None if planet_data.tidal_locked else planet_data.rotation_period,
        seasonal_period=None,   # fill later if orbital period is added
    )
