from src.input_data import set_env_parser
from src.input_data.input_configs import EnvironmentInputConfig
from src.input_data.planet_input import PlanetInputBuilder


def test_earth_planet_preset_builds():
    params = PlanetInputBuilder().build_planetary_params("Earth", None)
    star = PlanetInputBuilder().build_stellar_params("Earth", None)

    assert params.name == "Earth"
    assert star.name == "Sun"


def test_explicit_volcanic_rock_environment_input_is_preserved():
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

    environment_config = EnvironmentInputConfig.from_raw(
        {
            "environment": "volcanic_rock",
            "environment_source": "explicit",
            "environment_data": env_data,
        }
    )
    parsed = set_env_parser.build_local_env_data(
        environment_config=environment_config,
        planet_model="Earth",
    )

    assert parsed == env_data
