from src.environment.volcanic_rocks.volcanic_rock_env import VolcanicRock
from src.environment.hydro_env.hydrothermal_env import HydrothermalVent

#
#   set up environment object
#

def set_simulation_environment(protocell_list, env_type, env_data):
    # env_type = hydro_vent | volcanic_rock
    if env_type == "hydro_vent":
        return HydrothermalVent(protocell_list, env_data)
    elif env_type == "volcanic_rock":
        return VolcanicRock(protocell_list, env_data)
    elif env_type == "Titan":
        return VolcanicRock(protocell_list, env_data)
    else:
        raise NotImplementedError("simulation environment")
