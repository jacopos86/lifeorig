from src.environment.solvent import SolventData
from src.input_data.environment_base_parser import LocalEnvironmentParser
from src.utilities.logging_module import log
from src.environment.external_drive_params import ExternalDriveParams
from src.environment.liquid_level_model import LiquidLevelParams

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
            "pore_radius": self._quantity(
                input_dict=data.get("pore_radius"),
                required_keys=("units", "value"),
                desc="pore radius"
            ),
            "pore_height": self._quantity(
                input_dict=data.get("pore_height"),
                required_keys=("units", "value"),
                desc="pore height"
            ),
            "distance_neigh_pores": self._quantity(
                input_dict=data.get("distance_neigh_pores"),
                required_keys=("units", "value"),
                desc="distance neigh. pores"
            ),
            "temperature": self._quantity(
                input_dict=data.get("temperature"),
                required_keys=("units", "value"),
                desc="temperature"
            ),
            "pressure": self._quantity(
                input_dict=data.get("pressure"),
                required_keys=("units", "value"),
                desc="pressure"
            ),
            "solvent_data": self._get_solvent_data(data),
            "liquid_level_params": self._get_liquid_level_params(data)
        }
    def _get_solvent_data(self, env_data: dict) -> SolventData:
        data = env_data.get("solvent_data", {})
        return SolventData(
            name=data.get("name", "H2O"),
            density=self._quantity(
                input_dict=data.get("density"),
                required_keys=("units", "value"),
                desc="solvent density"
            ),
            dynamic_viscosity=self._quantity(
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
    def _get_liquid_level_params(self, env_data):
        data = env_data.get("solvent_data", {}).get("liquid_level_params", {})
        return LiquidLevelParams(
            model_type=data.get("type", "constant"),
            base_level=self._quantity(
                input_dict=data.get("base_level"),
                required_keys=("units", "value"),
                desc="liquid base level"
            )
        )
    def _parse_external_drive_params(self, data: dict, desc: str) -> ExternalDriveParams:
        return None
    def set_external_env_forces(self):
        return None