import math
from src.common.units import Q_
from src.environment.external_drive_params import TimeDependentField, LiquidLevelParams
from src.chemical_env.chemical_environment import ChemEnvResult
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.utilities.logging_module import log

#
#   Liquid level models
#

class ConstantLiquidLevel(TimeDependentField):
    def __init__(self, params: LiquidLevelParams):
        super().__init__(params)
    def value(self, t: Q_) -> Q_:
        return self.params.base_level


class SinusoidalLiquidLevel(TimeDependentField):
    def __init__(self, params: LiquidLevelParams):
        super().__init__(params)
    def value(self, t: Q_) -> Q_:
        return self.params.base_level + self.params.amplitude * math.sin(
            2.0 * math.pi * t / self.params.period + self.params.phase
        )


class PiecewiseLiquidLevel(TimeDependentField):
    def __init__(self, params: LiquidLevelParams):
        super().__init__(params)
    def value(self, t: Q_) -> Q_:
        for idx, switch_time in enumerate(self.params.switch_times):
            if t < switch_time:
                return self.params.levels[idx]
        return self.params.levels[-1]

#
#   build liquid level field
#

def build_liquid_level_field(params: LiquidLevelParams) -> TimeDependentField:
    model_type = params.model_type.lower()
    if model_type == "constant":
        return ConstantLiquidLevel(params)
    if model_type == "sinusoidal":
        return SinusoidalLiquidLevel(params)
    if model_type == "piecewise":
        return PiecewiseLiquidLevel(params)
    raise log.error(f"Unknown liquid level model_type: {params.model_type}")

#
#   derive chemical water factor
#

def derive_liquid_base_factor(
        env_model: str,
        env_data: dict,
        planet_data: PlanetaryEnvironmentParams,
        chemical_env: ChemEnvResult,
        vesc: Q_,
        min_factor: float = 0.0,
        max_factor: float = 2.0
    ) -> float:
    # check env model
    if env_model == "volcanic_rock":
        return derive_volcanic_rock_water_base_factor(
            chemical_env=chemical_env,
            vesc=vesc,
            min_factor=min_factor,
            max_factor=max_factor
        )
    if env_model == "hydro_vent":
        return None
    if env_model == "Titan":
        return derive_Titan_solvent_base_factor(
            planet_data=planet_data,
            env_data=env_data
        )
    log.error(f"Unknown environment model: {env_model}")

#
#    liquid amplitude factor
#

def derive_liquid_amplitude_factor(
        env_model: str,
        planet_data: PlanetaryEnvironmentParams,
        water_base_factor: float,
        min_factor: float = 0.0,
        max_factor: float = 2.0
    ) -> float:
    if env_model == "hydro_vent":
        return None
    if env_model == "volcanic_rock":
        return derive_volcanic_rock_water_amplitude_factor(
            planet_data=planet_data,
            water_base_factor=water_base_factor,
            min_factor=min_factor,
            max_factor=max_factor
        )
    log.error(f"Unknown environment model: {env_model}")