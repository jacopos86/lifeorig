from src.input_data.parse_volcanic_rock_env import VolcanicRockLocalEnvironmentParser
from src.input_data.parse_titan_env import TitanLocalEnvironmentParser
from src.utilities.logging_module import log

#
#    build local environment
#

def build_local_env_data(input_parser):
    env_model = input_parser.environment_config.env_type
    planet_model = input_parser._data.get("planet_model")
    if input_parser.environment_config.is_explicit:
        return input_parser.environment_config.data
    if planet_model == "Titan":
        return TitanLocalEnvironmentParser(
            input_parser=input_parser,
            env_data=input_parser._data.get("environment_data")
        )._build_local_env()
    if env_model == "volcanic_rock":
        return VolcanicRockLocalEnvironmentParser(
            input_parser=input_parser,
            env_data=input_parser._data.get("environment_data")
        )._build_local_env()
    if env_model == "hydro_vent":
        # TODO
        return None
    if env_model == "impact_crater":
        # TODO
        return None
    log.error(f"env model not recognized: {env_model}")