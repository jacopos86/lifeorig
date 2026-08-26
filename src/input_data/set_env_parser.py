from src.input_data.parse_volcanic_rock_env import VolcanicRockLocalEnvironmentParser
from src.input_data.parse_titan_env import TitanLocalEnvironmentParser
from src.utilities.logging_module import log

#
#    build local environment
#

def build_local_env_data(
    environment_config,
    planet_model,
    working_dir=None,
):
    env_model = environment_config.env_type
    env_data = environment_config.data
    if environment_config.is_explicit:
        return env_data
    if planet_model == "Titan":
        return TitanLocalEnvironmentParser(
            env_data=env_data,
            working_dir=working_dir,
        )._build_local_env()
    if env_model == "volcanic_rock":
        return VolcanicRockLocalEnvironmentParser(
            env_data=env_data,
            working_dir=working_dir,
        )._build_local_env()
    if env_model == "hydro_vent":
        # TODO
        return None
    if env_model == "impact_crater":
        # TODO
        return None
    log.error(f"env model not recognized: {env_model}")