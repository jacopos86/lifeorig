from src.environment.solvent import SolventData
from src.input_data.environment_base_parser import LocalEnvironmentParser
from src.utilities.logging_module import log

#
#   Volcanic Rock Local environment parser
#

class VolcanicRockLocalEnvironmentParser(LocalEnvironmentParser):
    def _build_local_env(self):
        data = self.env_data
        if data is None:
            log.error("Missing 'environment_data' in input")
        return {
            "num_pores": data.get("number_pores"),
            "pore_radius": self._parse_quantity_from_dict(
                input_dict=data.get("pore_radius"),
                required_keys=("units", "value"),
                desc="pore radius"
            ),
            "pore_height": self._parse_quantity_from_dict(
                input_dict=data.get("pore_height"),
                required_keys=("units", "value"),
                desc="pore height"
            ),
            "distance_neigh_pores": self._parse_quantity_from_dict(
                input_dict=data.get("distance_neigh_pores"),
                required_keys=("units", "value"),
                desc="distance neigh. pores"
            ),
            "temperature": self._parse_quantity_from_dict(
                input_dict=data.get("temperature"),
                required_keys=("units", "value"),
                desc="temperature"
            ),
            "pressure": self._parse_quantity_from_dict(
                input_dict=data.get("pressure"),
                required_keys=("units", "value"),
                desc="pressure"
            ),
            "solvent_data": self._get_solvent_data(data),
        }

    def _get_solvent_data(self, env_data: dict) -> SolventData:
        data = env_data.get("solvent_data", {})
        return SolventData(
            name=data.get("name", "H2O"),
            liquid_level_params=self._get_liquid_level_params(env_data),
            density=self._parse_quantity_from_dict(
                input_dict=data.get("density"),
                required_keys=("units", "value"),
                desc="solvent density"
            ),
            dynamic_viscosity=self._parse_quantity_from_dict(
                input_dict=data.get("dynamic_viscosity"),
                required_keys=("units", "value"),
                desc="solvent dynamic viscosity"
            ),
            dielectric_constant=(
                None if data.get("dielectric_constant") is None
                else float(data.get("dielectric_constant"))
            ),
            diffusion_scale=float(data.get("diffusion_scale", 1.0)),
            polarity=(
                None if data.get("polarity") is None
                else float(data.get("polarity"))
            ),
        )

    def set_external_env_forces(self):
        return None
