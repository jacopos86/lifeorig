from src.environment.volcanic_rocks.volcanic_rock_env import VolcanicRock
from src.environment.hydro_env.hydrothermal_env import HydrothermalVent
from src.environment.surf_pond.surface_pond_env import SurfacePondEnvironment

#
#   set up environment object
#

def setup_environment(input_params):
    # TITAN dispatcher
    if input_params.planetary_data.name == "Titan":
        return setup_titan_environment(input_params)
    # volcanic rock model case
    if input_params.env_model == "volcanic_rock":
        return setup_volcanic_rock_environment(input_params)
    # surface pond case
    if input_params.env_model == "surface_pond":
        return setup_surface_pond_environment(input_params)
    raise NotImplementedError("Environment not recognized")

#
#    set different environments
#

def setup_titan_environment(input_params):
    # planetary data
    planet = input_params.planetary_data
    # Star
    star = input_params.stellar_data
    # Local env properties
    local_env_data = input_params.local_env_data
    env_type = local_env_data["local_environment"]
    # set up local environment
    local_env = set_simulation_environment(
        env_type,
        local_env_data,
        star,
        planet,
        time_grid=input_params.time_grid,
        working_dir=input_params.working_dir,
    )
    exit()

    geometry = local["pool_geometry"]
    solvent = local["solvent_data"]
    external_forces = local["external_forces"]

    liquid_level_model = setup_titan_liquid_level(solvent, external_forces)
    radiation_model = setup_titan_radiation(planet, star, local)

    return {
        "type": "Titan_surface_pond",
        "planet": planet,
        "star": star,
        "geometry": geometry,
        "solvent": solvent,
        "external_forces": external_forces,
        "liquid_level_model": liquid_level_model,
        "radiation_model": radiation_model,
    }

#
#   call simulation environment
#

def set_simulation_environment(env_type, env_data, star_data, planet_data, time_grid=None, working_dir=None):
    # env_type = hydro_vent | volcanic_rock | surface pond
    if env_type == "hydro_vent":
        return HydrothermalVent(env_data)
    elif env_type == "surface_pond" or env_type == "methane_pool_column":
        return SurfacePondEnvironment(planet_data, star_data, env_data, time_grid=time_grid, working_dir=working_dir)
    elif env_type == "volcanic_rock":
        return VolcanicRock(env_data)
    else:
        raise NotImplementedError("simulation environment")
