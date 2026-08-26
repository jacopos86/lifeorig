from src.common.units import Q_
from src.environment.liquid_level_model import LiquidLevelParams
from src.environment.solvent import SolventData
from src.input_data.environment_base_parser import LocalEnvironmentParser
from src.utilities.logging_module import log
from src.utilities.plot_titan_pool_profile import plot_titan_surface_area_profile
from src.environment.external_drive_params import ExternalDriveParams

#
#   Titan environment parser
#

class TitanLocalEnvironmentParser(LocalEnvironmentParser):
    def _build_local_env(self):
        data = self.env_data
        pool_geometry = None
        if data is not None:
            pool_geometry = self._get_pool_geometry(data)
        #
        solvent_data = data.get("solvent_data", {})
        liquid_level_data = data.get("liquid_level_params", {})
        base_level = self._quantity(
            input_dict=liquid_level_data.get("base_level"),
            required_keys=("units", "value"),
            desc="methane liquid base level"
        )
        initial_level = self._quantity(
            input_dict=liquid_level_data.get("initial_level"),
            required_keys=("units", "value"),
            desc="methane liquid initial level"
        )
        methane_level_params = LiquidLevelParams(
            base_level=base_level,
            min_level=Q_(0.0, pool_geometry["height"].units),
            max_level=pool_geometry["height"],
            initial_level=initial_level if initial_level is not None else base_level,
        )
        # set evaporation levels + rainfall fluxes
        parsed_forces = self.set_external_env_forces()
        return {
            "local_environment": "methane_pool_column",
            "pool_geometry": pool_geometry,
            "temperature": Q_(94.0, "kelvin"),
            "pressure": Q_(1.45, "bar"),
            "solvent_data": SolventData(
                name=solvent_data.get("name", "CH4 based mixture"),
                composition=solvent_data.get("composition", {}),
                density=Q_(450.0, "kg / m^3"),
                dynamic_viscosity=Q_(1.8e-4, "Pa * s"),
                dielectric_constant=1.7,
                diffusion_scale=1.0,
                polarity=0.0,
            ),
            "liquid_level_params": methane_level_params,
            "external_forces": parsed_forces
        }
    # get geometry
    def _get_pool_geometry(self, env_data: dict):
        geometry_data = env_data.get("pool_geometry")
        if geometry_data is None:
            return None
        height = self._quantity(
            input_dict=geometry_data.get("height"),
            required_keys=("units", "value"),
            desc="Titan pool height"
        )
        surface_area_z0 = self._quantity(
            input_dict=geometry_data.get("surface_area_z0"),
            required_keys=("units", "value"),
            desc="Titan pool bottom surface area"
        )
        # R_A / S_A ratio
        reactive_area_surface_area_ratio = geometry_data.get("reactive_area_to_surface_area_ratio")
        if reactive_area_surface_area_ratio is None:
            log.error("Missing pool_geometry.reactive_area_to_surface_area_ratio")
        # S_A(z)
        surface_area_profile = self._get_surface_area_profile(geometry_data)
        surface_area_profile_plot = plot_titan_surface_area_profile(
            output_dir=self.working_dir,
            height=height,
            surface_area=surface_area_z0,
            profile_data=surface_area_profile,
            reactive_area_to_surface_area_ratio=float(reactive_area_surface_area_ratio),
        )
        return {
            "height": height,
            "n_grid_cells": int(geometry_data.get("n_grid_cells", 100)),
            "surface_area_z0": surface_area_z0,
            "surface_area_profile": surface_area_profile,
            "reactive_area_to_surface_area_ratio": float(reactive_area_surface_area_ratio),
            "surface_area_profile_plot": surface_area_profile_plot
        }
    # get profile S_A(z)
    def _get_surface_area_profile(self, geometry_data: dict):
        profile_data = geometry_data.get("surface_area_profile", {"type": "uniform"})
        profile_type = profile_data.get("type", "uniform")
        profile = {"type": profile_type}
        if profile_type == "bottom_weighted":
            profile["decay_length"] = self._quantity(
                input_dict=profile_data.get("decay_length"),
                required_keys=("units", "value"),
                desc="Titan surface area profile decay length"
            )
        elif profile_type != "uniform":
            log.error(f"Unknown Titan surface area profile type: {profile_type}")
        return profile
    # get solvent data
    def _get_solvent_data(self, env_data: dict) -> SolventData:
        data = env_data.get("solvent_data", {})
        return SolventData(
            name=data.get("name"),
            composition=data.get("composition", {}),
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
    # parse external drive
    def _parse_external_drive_params(self, data: dict, desc: str) -> ExternalDriveParams:
        model_type = data.get("model_type", data.get("type"))
        # base level
        base_level = self._quantity(
            input_dict=data.get("base_level"),
            required_keys=("units", "value"),
            desc=f"{desc} base level"
        )
        # amplitude
        amplitude = self._quantity(
            input_dict=data.get("amplitude"),
            required_keys=("units", "value"),
            desc=f"{desc} amplitude"
        )
        # period
        period = self._quantity(
            input_dict=data.get("period"),
            required_keys=("units", "value"),
            desc=f"{desc} period"
        )
        # switch times
        switch_times = self._quantity(
            input_dict=data.get("switch_times"),
            required_keys=("units", "value"),
            desc=f"{desc} switch times"
        )
        # levels
        levels = self._quantity(
            input_dict=data.get("levels"),
            required_keys=("units", "value"),
            desc=f"{desc} levels"
        )
        return ExternalDriveParams(
            model_type=model_type,
            base_level=base_level,
            amplitude=amplitude,
            period=period,
            phase=float(data.get("phase", 0.0)),
            switch_times=switch_times,
            levels=levels,
        )
    #
    # external forces set up driver
    def set_external_env_forces(self):
        data = self.env_data.get("external_forces", {})
        # set rainfall / evaporation params
        rainfall_params = None
        evaporation_params = None
        if "rainfall" in data:
            rainfall_params = self._parse_external_drive_params(
                data=data["rainfall"],
                desc="rainfall"
            )
        if "evaporation" in data:
            evaporation_params = self._parse_external_drive_params(
                data=data["evaporation"],
                desc="evaporation"
            )
        return {
            "rainfall": rainfall_params,
            "evaporation": evaporation_params
        }
