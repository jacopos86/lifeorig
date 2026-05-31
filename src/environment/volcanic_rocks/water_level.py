import numpy as np
from src.common.units import Q_
from src.utilities.useful_func import _clamp
from src.chemical_env.chemical_environment import ChemEnvResult
from src.planet_params.planetary_params import PlanetaryEnvironmentParams

#
#   water base factor
#

def derive_volcanic_rock_water_base_factor(
        chemical_env: ChemEnvResult,
        vesc: Q_,
        min_factor: float = 0.0,
        max_factor: float = 2.0
    ) -> float:
    # atmospheric H2O factor
    chem_water_factor = derive_atmospheric_h2o_water_factor(
        chemical_env=chemical_env
    )
    retention_factor = derive_retention_factor_from_escape_velocity(
        vesc=vesc
    )
    return _clamp(
        chem_water_factor * retention_factor,
        min_factor,
        max_factor
    )

def derive_atmospheric_h2o_water_factor(
        chemical_env: ChemEnvResult,
        species: str = "H2O",
        ref_h2o_partial_pressure: Q_ = Q_(0.01, "bar"),
        log_scale: float = 0.3,
        min_factor: float = 0.0,
        max_factor: float = 2.0
    ) -> float:
    if chemical_env is None:
        log.error("chemical environment cannot be None")
    h2o_partial_pressure = chemical_env.get_partial_pressure(species, default=0.0)
    ratio = (
        h2o_partial_pressure.to("bar").magnitude
        / ref_h2o_partial_pressure.to("bar").magnitude
    )
    if ratio <= 0.0:
        return 0.0
    factor = 1.0 + log_scale * np.log10(ratio)
    return _clamp(factor, min_factor, max_factor)

def derive_retention_factor_from_escape_velocity(
        vesc: Q_,
        ref_vesc: Q_ = Q_(11185.98, "m / s"),
        min_factor: float = 0.2,
        max_factor: float = 1.5
    ) -> float:
    ratio = (
        vesc.to("m / s").magnitude
         / ref_vesc.to("m / s").magnitude
    )
    return _clamp(ratio, min_factor, max_factor)

#
#   derive water amplitude factors
#

def derive_volcanic_rock_water_amplitude_factor(
        planet_data: PlanetaryEnvironmentParams,
        water_base_factor: float,
        min_factor: float,
        max_factor: float
    ) -> float:
    # eccentricity
    eccentricity_factor = 1.0 + 2.0 * planet_data.eccentricity
    # contrast day / night
    contrast_factor = 1.0 + planet_data.day_night_contrast
    if planet_data.tidal_locked:
        contrast_factor *= 2.0
    water_avail_factor = min(1.0, water_base_factor)
    # amplitude factor
    amplitude_factor = (
        eccentricity_factor
        * contrast_factor
        * water_avail_factor
    )
    return _clamp(amplitude_factor, min_factor, max_factor)