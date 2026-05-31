import os
from src.utilities.logging_module import log
from src.input_data.core import AbstractInput
from src.planet_params.planetary_params import PlanetaryEnvironmentParams
from src.planet_params.earth_params import get_earth_planetary_params, get_earth_stellar_params
from src.planet_params.europa_params import get_europa_planetary_params, get_europa_stellar_params
from src.planet_params.mars_params import get_mars_planetary_params, get_mars_stellar_params
from src.planet_params.titan_params import get_titan_planetary_params, get_titan_stellar_params
from src.planet_params.venus_params import get_venus_planetary_params, get_venus_stellar_params
from src.stellar_params.stellar_data import StellarParams
from src.environment.external_drive_params import LiquidLevelParams
from src.environment.solvent import SolventData
from src.environment.set_planet_environment import derive_planet_env_data
from src.chemical_env.chemical_environment import set_ChemEnvInput, SimulateChemEnv
from src.common.units import Q_
from src.environment.pool_spatial_profile import get_pool_spatial_profile_function
from src.utilities.plot_titan_pool_profile import plot_titan_pool_profile

#
#  parameters class

class parameters_class(AbstractInput):
    ''' parameters class '''
    _ALLOWED_METABOLITE_TYPES = {"binary", "multi"}
    _ALLOWED_METABOLITE_DISTR = {"uniform", "length_decay"}
    _ALLOWED_ENV_MODELS = {"volcanic_rock", "hydro_vent", "impact_crater"}
    _ALLOWED_PLANET_MODELS = {"Earth", "Europa", "Mars", "Titan", "Venus", "custom", None}
    _REQUIRED_ENV_VOLCROCK_KEYS = {
        "num_pores", 
        "pore_radius", 
        "pore_height", 
        "distance_neigh_pores", 
        "solvent_data"
    }
    def __init__(self):
        # work dir
        self.working_dir = None
        # n. protocell (initial)
        self.QSP_size = None
        # protocell info
        self.protocell_info = None
        # catalyst set parameters
        self.catalyst_set_params = None
        # rates distribution parameters
        self.rates_params = None
        # molecules data parameters
        self.metabolites_params = None
        # env. data
        self.env_model = None
        self.local_env_data = None
        # planetary data
        self.planetary_data = {}
        self.chem_env_params = None
        # stellar data
        self.stellar_data = None
    def _parse_data(self):
        # read input parameters
        if "working_dir" in self._data:
            self.working_dir = self._data["working_dir"]
            isExist = os.path.exists(self.working_dir)
            if not isExist:
                os.mkdir(self.working_dir)
        # planetary parameters
        planet_params = self.build_planetary_params(
            planet_model=self._data.get("planet_model"),
            planetary_data=self._data.get("planetary_data")
        )
        if planet_params is not None:
            self.planetary_data["basic_info"] = planet_params
            self.planetary_data["basic_info"].log_summary()
        # stellar parameters
        self.stellar_data = self.build_stellar_params(
            planet_model=self._data.get("planet_model"),
            stellar_data=self._data.get("stellar_data")
        )
        self.stellar_data.log_summary()
        # local environment data
        self.set_local_environment()
        # num. individuals in QSP to average
        if "QSP_size" in self._data:
            self.QSP_size = self._data["QSP_size"]
        # protocell parameters
        if "protocell_data" in self._data:
            self.protocell_info = {
                "n_shells": self._data["protocell_data"].get("n_shells"),
                "baseline_radius": self._parse_quantity_from_dict(
                    input_dict= self._data["protocell_data"].get("baseline_radius"), 
                    required_keys=("units", "value"), 
                    desc="protocell baseline radius") 
            }
        # metabolites data
        if "metabolites_data" in self._data:
            self.metabolites_params = self._data["metabolites_data"]
        # catalysts set size
        if "catalyst_set" in self._data:
            self.catalyst_set_params = self._data["catalyst_set"]
        #
        # mutation parameters
        #
        if "distribution_rates" in self._data:
            self.rates_params = self._data["distribution_rates"]
        # time variables
        # size sample space
        if "evol_params" in self._data:
            self.evol_size = self._data["evol_params"]
    # set local environment data and derived planetary forcing
    def set_local_environment(self):
        self.env_model = self._data.get("environment")
        self.local_env_data = self._build_local_env_data()
        if self.local_env_data is None:
            return
        # set input Chem. Env.
        chem_input = set_ChemEnvInput(
            planet_data=self.planetary_data,
            env_data=self.local_env_data
        )
        # simulate chem. structure environment
        self.chem_env_params = SimulateChemEnv(
            chem_input=chem_input,
            planet_data=self.planetary_data["basic_info"],
            stellar_data=self.stellar_data,
            atmosphere_data=self.planetary_data["basic_info"].atmosphere,
            output_dir=self.working_dir
        ).run()
        exit()
        # set derived planetary data
        derived_planet_data = derive_planet_env_data(
            self.env_model,
            self.local_env_data,
            self.planetary_data["basic_info"],
            self.chem_env_params
        )
        self.planetary_data["derived_params"] = derived_planet_data
        self._update_env_data_based_on_planet_data()
    # update environment data from derived planetary data
    def _update_env_data_based_on_planet_data(self):
        derived_params = self.planetary_data.get("derived_params")
        if derived_params is None:
            return
        self.local_env_data["gravity"] = derived_params.gravity
        solvent_data = self.local_env_data.get("solvent_data")
        liquid_level_params = solvent_data.liquid_level_params if solvent_data is not None else None
        if liquid_level_params is not None:
            if liquid_level_params.base_level is None:
                liquid_level_params.base_level = derived_params.water_base_factor
            if liquid_level_params.amplitude is None:
                liquid_level_params.amplitude = derived_params.water_amplitude_factor
            if liquid_level_params.period is None:
                liquid_level_params.period = derived_params.day_night_period
    # build local environment parameters
    def _build_local_env_data(self):
        env_model = self._data.get("environment")
        planet_model = self._data.get("planet_model")
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
            input_dict= data.get("base_level"), 
            required_keys=("units", "value"), 
            desc="liquid base level"
        )
        amplitude = self._parse_quantity_from_dict(
            input_dict= data.get("amplitude"), 
            required_keys=("units", "value"), 
            desc="liquid amplitude level"
        )
        period = self._parse_quantity_from_dict(
            input_dict= data.get("period"), 
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
    #
    #    build planetary parameters from input data
    #
    def build_planetary_params(
        self,
        planet_model: str,
        planetary_data: dict | None,
    ) -> PlanetaryEnvironmentParams:
        if planet_model == "Earth":
            return get_earth_planetary_params()
        if planet_model == "Europa":
            return get_europa_planetary_params()
        if planet_model == "Mars":
            return get_mars_planetary_params()
        if planet_model == "Titan":
            return get_titan_planetary_params()
        if planet_model == "Venus":
            return get_venus_planetary_params()
        if planetary_data is None:
            return None
        # return planetary data
        return PlanetaryEnvironmentParams(
            name=planetary_data.get("name", planet_model or "custom"),
            planet_radius=self._parse_quantity_from_dict(
                input_dict=planetary_data.get("planet_radius"),
                required_keys=("units", "value"), 
                desc="planet radius"
            ),
            planet_mass=self._parse_quantity_from_dict(
                input_dict=planetary_data.get("planet_mass"),
                required_keys=("units", "value"), 
                desc="planet mass"
            ),
            orbital_distance=self._parse_quantity_from_dict(
                input_dict=planetary_data.get("orbital_distance"),
                required_keys=("units", "value"), 
                desc="orbital distance"
            ),
            rotation_period=self._parse_quantity_from_dict(
                input_dict=planetary_data.get("rotation_period"),
                required_keys=("units", "value"), 
                desc="planet rotation period",
            ),
            obliquity=float(planetary_data.get("obliquity", 0.0)),
            eccentricity=float(planetary_data.get("eccentricity", 0.0)),
            tidal_locked=bool(planetary_data.get("tidal_locked", False)),
            day_night_contrast=float(planetary_data.get("day_night_contrast", 0.0)),
            chemical_env=planetary_data.get("exo_chemistry"),
            atmosphere=self._build_atmosphere_params(planetary_data.get("atmosphere"))
        )
    # build atmospheric data
    def _build_atmosphere_params(self, atmosphere_data: dict | None) -> dict | None:
        if atmosphere_data is None:
            return None
        return {
            "n_layers": atmosphere_data.get("n_layers"),
            "z_max": self._parse_quantity_from_dict(
                input_dict=atmosphere_data.get("z_max"),
                required_keys=("units", "value"),
                desc="atmosphere z_max"
            ),
            "top_pressure": self._parse_quantity_from_dict(
                input_dict=atmosphere_data.get("top_pressure"),
                required_keys=("units", "value"),
                desc="atmosphere top_pressure"
            ),
        }
    #
    #    build stellar parameters from input data
    #
    def build_stellar_params(
        self,
        planet_model: str,
        stellar_data: dict | None,
    ) -> StellarParams | None:
        if planet_model == "Earth":
            return get_earth_stellar_params()
        if planet_model == "Europa":
            return get_europa_stellar_params()
        if planet_model == "Mars":
            return get_mars_stellar_params()
        if planet_model == "Titan":
            return get_titan_stellar_params()
        if planet_model == "Venus":
            return get_venus_stellar_params()
        if stellar_data is None:
            return None
        return StellarParams(
            name=stellar_data.get("name", "star"),
            spectral_class=stellar_data.get("spectral_class"),
            effective_temperature=self._parse_quantity_from_dict(
                input_dict=stellar_data.get("star_temperature"),
                required_keys=("units", "value"),
                desc="star temperature"
            ),
            radius=self._parse_quantity_from_dict(
                input_dict=stellar_data.get("star_radius"),
                required_keys=("units", "value"),
                desc="star radius"
            ),
            mass=self._parse_quantity_from_dict(
                input_dict=stellar_data.get("star_mass"),
                required_keys=("units", "value"),
                desc="star mass"
            ),
            luminosity=self._parse_quantity_from_dict(
                input_dict=stellar_data.get("luminosity"),
                required_keys=("units", "value"),
                desc="star luminosity"
            ),
        )
    #
    #    validation section
    #
    def _validate(self):
        # optional: check metabolites_parameters
        required_keys = ["type", "pol_strng_maxsize", "metabolites_distr_type", "initial_population_molecules"]
        missing = [k for k in required_keys 
           if k not in self.metabolites_params or self.metabolites_params[k] is None]
        if missing:
            log.error(f"Missing keys in metabolites_params: {missing}")
        # check type
        if self.metabolites_params.get("type") not in self._ALLOWED_METABOLITE_TYPES:
            log.error(
                f"Invalid type '{self.metabolites_params.get('type')}'. "
                f"Valid options: {sorted(self._ALLOWED_METABOLITE_TYPES)}"
        )
        # check distr. model
        if self.metabolites_params.get("metabolites_distr_type") not in self._ALLOWED_METABOLITE_DISTR:
            log.error(
                f"Invalid distr. type '{self.metabolites_params.get('metabolites_distr_type')}'. "
                f"Valid options: {sorted(self._ALLOWED_METABOLITE_DISTR)}"
            )
        # check environment model
        if self.env_model not in self._ALLOWED_ENV_MODELS:
            log.error(
                f"Invalid env. type '{self.env_model}'. "
                f"Valid options: {sorted(self._ALLOWED_ENV_MODELS)}"
            )
        self._validate_planet_data()
        self._validate_environment_data()
    # validate planet input
    def _validate_planet_data(self):
        planet_model = self._data.get("planet_model")
        planetary_data = self._data.get("planetary_data")
        if planet_model not in self._ALLOWED_PLANET_MODELS:
            log.error(
                f"Invalid planet_model '{planet_model}'. "
                f"Valid options: {sorted(model for model in self._ALLOWED_PLANET_MODELS if model is not None)}"
            )
        if planet_model in {"custom", None} and planetary_data is None:
            log.error("planetary_data is required when planet_model is custom or missing")
    # validate environment data
    def _validate_environment_data(self):
        """
        Validate environment input dictionary.
        Raises ValueError if something is wrong.
        """
        if self.env_model == "volcanic_rock":
            # check keys
            missing = self._REQUIRED_ENV_VOLCROCK_KEYS - set(self.env_data.keys())
            if missing:
                log.error(f"Missing environment keys: {missing}")
            # number of pores
            n_pores = self.env_data["num_pores"]
            if not isinstance(n_pores, int) or n_pores <= 0:
                log.error("number_pores must be a positive integer")
            # pore geometry
            pore_radius = self.env_data["pore_radius"]
            pore_height = self.env_data["pore_height"]
            if pore_radius <= 0:
                log.error("pore_radius must be > 0")
            if pore_height <= 0:
                log.error("pore_height must be > 0")
            # distance between closest neighbors
            min_distance = self.env_data["distance_neigh_pores"]
            if min_distance <= 0:
                log.error("min_distance must be > 0")
        # ----------------------------
        # temperature / pressure
        # ----------------------------
        T = self.env_data["temperature"]
        P = self.env_data["pressure"]
        if T <= 0:
            log.error("temperature must be > 0")
        if P <= 0:
            log.error("pressure must be > 0")
        self._validate_liquid_level_params()
    # validate liquid level parameters
    def _validate_liquid_level_params(self):
        solvent_data = self.env_data.get("solvent_data")
        params = solvent_data.liquid_level_params if solvent_data is not None else None
        if params is None:
            log.error("Missing solvent_data.liquid_level_params in environment data")
        if params.model_type not in {"constant", "sinusoidal", "piecewise"}:
            log.error(f"Unknown liquid level model_type: {params.model_type}")
        if params.base_level is None:
            log.error("liquid_level_params: 'base_level' is required")
        if params.period is not None and params.period <= 0:
            log.error("liquid_level_params: 'period' must be > 0 when provided")

p = parameters_class()
p.sep = "*"*94
