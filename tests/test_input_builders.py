from src.input_data.environment_input import EnvironmentInputBuilder
from src.input_data.core import AbstractInput
from src.input_data.input_configs import EnvironmentInputConfig
from src.input_data.planet_input import PlanetInputBuilder


class _EnvironmentBuilderHarness(AbstractInput, EnvironmentInputBuilder, PlanetInputBuilder):
    def __init__(self, env_data):
        self._data = {
            "planet_model": "Earth",
            "environment": "volcanic_rock",
            "environment_data": env_data,
        }
        self.environment_config = EnvironmentInputConfig.from_raw(self._data)

    def _parse_data(self):
        raise NotImplementedError

    def _validate(self):
        raise NotImplementedError


def test_earth_planet_preset_builds():
    params = PlanetInputBuilder().build_planetary_params("Earth", None)
    star = PlanetInputBuilder().build_stellar_params("Earth", None)

    assert params.name == "Earth"
    assert star.name == "Sun"


def test_volcanic_rock_environment_input_parses_quantities():
    env_data = {
        "number_pores": 10,
        "pore_radius": {"value": 1.0, "units": "millimeter"},
        "pore_height": {"value": 10.0, "units": "millimeter"},
        "distance_neigh_pores": {"value": 5.0, "units": "millimeter"},
        "temperature": {"value": 298.15, "units": "kelvin"},
        "pressure": {"value": 1.0, "units": "bar"},
        "solvent_data": {
            "name": "H2O",
            "liquid_level_params": {
                "type": "constant",
                "base_level": {"value": 1.0, "units": "millimeter"},
            },
            "density": {"value": 997.0, "units": "kg / m^3"},
            "dynamic_viscosity": {"value": 8.9e-4, "units": "Pa * s"},
            "dielectric_constant": 78.4,
            "diffusion_scale": 1.0,
            "polarity": 1.0,
        },
    }

    parsed = _EnvironmentBuilderHarness(env_data)._build_local_env_data()

    assert parsed["num_pores"] == 10
    assert parsed["pore_radius"].magnitude == 1.0
    assert parsed["solvent_data"].name == "H2O"
    assert parsed["solvent_data"].liquid_level_params.model_type == "constant"
