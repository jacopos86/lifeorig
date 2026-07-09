from src.common.units import Q_
from src.environment.external_drive_params import LiquidLevelParams
from src.environment.solvent import SolventData
from src.input_data.environment_base_parser import LocalEnvironmentParser
from src.utilities.logging_module import log
from src.utilities.plot_titan_pool_profile import plot_titan_surface_area_profile

#
#   Titan environment parser
#

class TitanLocalEnvironmentParser(LocalEnvironmentParser):
    def _build_local_env(self):
        data = self.env_data
        pool_geometry = None
        if data is not None:
            pool_geometry = self._get_pool_geometry(data)
        # methane solvent base level
        exit()
        methane_level_params = LiquidLevelParams(
            model_type="sinusoidal",
            base_level=Q_(5.0, "millimeter"),
            amplitude=Q_(2.0, "millimeter"),
            period=Q_(15.945, "day"),
            phase=0.0,
        )

        self.set_external_env_forces()

        return {
            "local_environment": "methane_pool",
            "pool_geometry": pool_geometry,
            "temperature": Q_(94.0, "kelvin"),
            "pressure": Q_(1.45, "bar"),
            "solvent_data": SolventData(
                name="CH4 based mixture",
                liquid_level_params=methane_level_params,
                density=Q_(450.0, "kg / m^3"),
                dynamic_viscosity=Q_(1.8e-4, "Pa * s"),
                dielectric_constant=1.7,
                diffusion_scale=1.0,
                polarity=0.0,
            ),
        }

    def _get_pool_geometry(self, env_data: dict):
        geometry_data = env_data.get("pool_geometry")
        if geometry_data is None:
            return None
        height = self._parse_quantity_from_dict(
            input_dict=geometry_data.get("height"),
            required_keys=("units", "value"),
            desc="Titan pool height"
        )
        surface_area = self._parse_quantity_from_dict(
            input_dict=geometry_data.get("surface_area"),
            required_keys=("units", "value"),
            desc="Titan pool surface area"
        )
        reactive_area_surface_area_ratio = geometry_data.get("reactive_area_surface_area_ratio")
        if reactive_area_surface_area_ratio is None:
            reactive_area_surface_area_ratio = geometry_data.get("reactive_area_to_surface_area_ratio")
        if reactive_area_surface_area_ratio is None:
            log.error("Missing pool_geometry.reactive_area_to_surface_area_ratio")
        surface_area_profile = self._get_surface_area_profile(geometry_data)
        surface_area_profile_plot = plot_titan_surface_area_profile(
            output_dir=self.working_dir,
            height=height,
            surface_area=surface_area,
            profile_data=surface_area_profile,
            reactive_area_to_surface_area_ratio=float(reactive_area_surface_area_ratio),
        )
        return {
            "height": height,
            "surface_area": surface_area,
            "surface_area_profile": surface_area_profile,
            "reactive_area_to_surface_area_ratio": float(reactive_area_surface_area_ratio),
            "reactive_area": surface_area * float(reactive_area_surface_area_ratio),
            "surface_area_profile_plot": surface_area_profile_plot,
        }

    def _get_surface_area_profile(self, geometry_data: dict):
        profile_data = geometry_data.get("surface_area_profile", {"type": "uniform"})
        profile_type = profile_data.get("type", "uniform")
        profile = {"type": profile_type}
        if profile_type == "bottom_weighted":
            profile["decay_length"] = self._parse_quantity_from_dict(
                input_dict=profile_data.get("decay_length"),
                required_keys=("units", "value"),
                desc="Titan surface area profile decay length"
            )
        elif profile_type != "uniform":
            log.error(f"Unknown Titan surface area profile type: {profile_type}")
        return profile

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
