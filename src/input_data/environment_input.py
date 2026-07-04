from src.common.units import Q_
from src.environment.external_drive_params import LiquidLevelParams
from src.environment.solvent import SolventData
from src.environment.surf_pond.pool_spatial_profile import get_pool_spatial_profile_function
from src.utilities.logging_module import log
from src.utilities.plot_titan_pool_profile import plot_titan_pool_profile

#
#   General Environment Builder
#

class EnvironmentInputBuilder:
    # build local environment parameters
    def _build_local_env_data(self):
        env_model = self.environment_config.env_type
        planet_model = self._data.get("planet_model")
        if self.environment_config.is_explicit:
            return self.environment_config.data
        if planet_model == "Titan":
            return self._build_titan_local_env_data()
        if env_model == "volcanic_rock":
            return self._build_volcanic_rock_env_data(self._data.get("environment_data"))
        if env_model == "hydro_vent":
            # TODO
            return None
        if env_model == "impact_crater":
            # TODO
            return None
        log.error(f"env model not recognized: {env_model}")

    # build volcanic rock local environment data
    def _build_volcanic_rock_env_data(self, data: dict | None):
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

    # build Titan local environment data
    def _build_titan_local_env_data(self):
        data = self._data.get("environment_data")
        pool_spatial_profile = None
        if data is not None:
            pool_spatial_profile = self._get_pool_spatial_profile(data)
            if pool_spatial_profile is not None:
                plot_titan_pool_profile(
                    output_dir=self.working_dir,
                    function_name=pool_spatial_profile["function_name"],
                    length=pool_spatial_profile["domain_length"],
                )
        # methane level
        exit()
        methane_level_params = LiquidLevelParams(
            model_type="sinusoidal",
            base_level=Q_(5.0, "millimeter"),
            amplitude=Q_(2.0, "millimeter"),
            period=Q_(15.945, "day"),
            phase=0.0,
        )
        return {
            "local_environment": "methane_pool",
            "pool_spatial_profile": pool_spatial_profile,
            "temperature": Q_(94.0, "kelvin"),
            "pressure": Q_(1.45, "bar"),
            "solvent_data": SolventData(
                name="CH4",
                liquid_level_params=methane_level_params,
                density=Q_(450.0, "kg / m^3"),
                dynamic_viscosity=Q_(1.8e-4, "Pa * s"),
                dielectric_constant=1.7,
                diffusion_scale=1.0,
                polarity=0.0,
            ),
        }

    # get pool spatial profile
    def _get_pool_spatial_profile(self, env_data: dict):
        profile_data = env_data.get("pool_spatial_profile")
        if profile_data is None:
            return None
        bottom_profile = profile_data.get("bottom_profile", {})
        function_name = bottom_profile.get("function")
        if function_name is None:
            log.error("Missing pool_spatial_profile.bottom_profile.function")
        length = self._parse_quantity_from_dict(
            input_dict=profile_data.get("domain_length"),
            required_keys=("units", "value"),
            desc="pool spatial profile domain length"
        )
        return {
            "type": profile_data.get("type"),
            "domain_length": length,
            "function_name": function_name,
            "function": get_pool_spatial_profile_function(function_name),
        }

    # get solvent data
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

    # get liquid level parameters
    def _get_liquid_level_params(self, env_data: dict) -> LiquidLevelParams:
        solvent_data = env_data.get("solvent_data", {})
        if "liquid_level_params" in solvent_data:
            data = solvent_data.get("liquid_level_params")
        elif "liquid_level_params" in env_data:
            data = env_data.get("liquid_level_params")
        else:
            data = env_data.get("water_level_params")
        if data is None:
            log.error("Missing 'liquid_level_params' in input")
        # --- required ---
        try:
            model_type = data.get("model_type", data.get("type")).lower()
        except KeyError as e:
            log.error(f"Missing key in liquid_level_params: {e}")
        except AttributeError:
            log.error("Missing key in liquid_level_params: model_type")
        # --- optional ---
        base_level = self._parse_quantity_from_dict(
            input_dict=data.get("base_level"),
            required_keys=("units", "value"),
            desc="liquid base level"
        )
        amplitude = self._parse_quantity_from_dict(
            input_dict=data.get("amplitude"),
            required_keys=("units", "value"),
            desc="liquid amplitude level"
        )
        period = self._parse_quantity_from_dict(
            input_dict=data.get("period"),
            required_keys=("units", "value"),
            desc="liquid oscillation period"
        )
        phase = float(data.get("phase", 0.0))
        #  piecewise model
        switch_times = self._parse_quantity_from_dict(
            input_dict=data.get("switch_times"),
            required_keys=("units", "value"),
            desc="liquid level switch times"
        )
        levels = self._parse_quantity_from_dict(
            input_dict=data.get("levels"),
            required_keys=("units", "value"),
            desc="piecewise water levels"
        )
        return LiquidLevelParams(
            model_type=model_type,
            base_level=base_level,
            amplitude=amplitude,
            period=period,
            phase=phase,
            switch_times=switch_times,
            levels=levels,
        )
