import os
from src.utilities.logging_module import log
from src.input_data.core import AbstractInput
from src.input_data.input_configs import (
    ChemicalNetworkInputConfig,
    EnvironmentInputConfig,
    EvolutionInputConfig,
    MoleculeInputConfig,
)
from src.input_data.set_env_parser import build_local_env_data
from src.input_data.planet_input import PlanetInputBuilder
from src.grids.temporal_grid import TimeGrid

#
#  parameters class

class parameters_class(AbstractInput, PlanetInputBuilder):
    ''' parameters class '''
    _ALLOWED_METABOLITE_TYPES = {"binary", "multi", "reference_file"}
    _ALLOWED_METABOLITE_DISTR = {"uniform", "length_decay"}
    _ALLOWED_ENV_MODELS = {"volcanic_rock", "hydro_vent", "surface_pond", "impact_crater"}
    _ALLOWED_ENV_SOURCES = {"planetary_solver", "explicit"}
    _ALLOWED_PLANET_MODELS = {"Earth", "Europa", "Mars", "Titan", "Venus", "custom", None}
    _REQUIRED_ENV_VOLCROCK_KEYS = {
        "num_pores", 
        "pore_radius", 
        "pore_height", 
        "distance_neigh_pores", 
        "solvent_data"
    }
    def __init__(self):
        self.sep = "*"*94
        # work dir
        self.working_dir = None
        # n. protocell (initial)
        self.QSP_size = None
        # protocell info
        self.protocell_info = None
        # catalyst set parameters
        self.catalyst_set_params = None
        # chemical network data
        self.chemical_network_data = None
        self.chemistry_config = None
        # rates distribution parameters
        self.rates_params = None
        # molecules data parameters
        self.metabolites_params = None
        self.molecule_config = None
        # env. data
        self.env_model = None
        self.environment_config = None
        self.local_env_data = None
        # planetary data
        self.planetary_data = None
        self.chem_env_params = None
        # stellar data
        self.stellar_data = None
        # evolution data
        self.evolution_config = None
        self.time_grid = None
    def _parse_data(self):
        # read input parameters
        self.environment_config = EnvironmentInputConfig.from_raw(self._data)
        self.chemistry_config = ChemicalNetworkInputConfig.from_raw(self._data)
        self.molecule_config = MoleculeInputConfig.from_raw(self._data)
        self.evolution_config = EvolutionInputConfig.from_raw(self._data)
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
            self.planetary_data = planet_params
            self.planetary_data.log_summary()
        # stellar parameters
        self.stellar_data = self.build_stellar_params(
            planet_model=self._data.get("planet_model"),
            stellar_data=self._data.get("stellar_data")
        )
        self.stellar_data.log_summary()
        # local environment data
        self.set_local_environment()
        # time grid
        if "time_grid" in self._data:
            self.time_grid = self._parse_time_grid(self._data["time_grid"])
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
            self.metabolites_params = self.molecule_config.data
        # catalysts set size
        if "catalyst_set" in self._data:
            self.catalyst_set_params = self._data["catalyst_set"]
        # chemical network
        if "chemical_network" in self._data:
            self.chemical_network_data = self.chemistry_config.data
        #
        # mutation parameters
        #
        if "distribution_rates" in self._data:
            self.rates_params = self._data["distribution_rates"]
        # time variables
        # size sample space
        if "evol_params" in self._data:
            self.evol_size = self.evolution_config.data
    # set local environment data
    def set_local_environment(self):
        self.env_model = self.environment_config.env_type
        self.local_env_data = build_local_env_data(
            environment_config=self.environment_config,
            planet_model=self._data.get("planet_model"),
            working_dir=self.working_dir,
        )
    #
    #    validation section
    #
    def _validate(self):
        # optional: check metabolites_parameters
        required_keys = ["type", "initial_population_molecules"]
        if self.metabolites_params.get("type") in {"binary", "multi"}:
            required_keys.extend(["pol_strng_maxsize", "metabolites_distr_type"])
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
        if (
            self.metabolites_params.get("type") in {"binary", "multi"}
            and self.metabolites_params.get("metabolites_distr_type") not in self._ALLOWED_METABOLITE_DISTR
        ):
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
        if self.environment_config.source not in self._ALLOWED_ENV_SOURCES:
            log.error(
                f"Invalid environment source '{self.environment_config.source}'. "
                f"Valid options: {sorted(self._ALLOWED_ENV_SOURCES)}"
            )
        self._validate_planet_data()
        self._validate_time_grid()
        self._validate_environment_data()
    # parse time grid
    def _parse_time_grid(self, data):
        start = self._parse_quantity_from_dict(data.get("start"), desc="time_grid.start")
        end = self._parse_quantity_from_dict(data.get("end"), desc="time_grid.end")
        dt = self._parse_quantity_from_dict(data.get("dt"), desc="time_grid.dt")
        nt = int(round(((end - start) / dt).to_base_units().magnitude))
        return TimeGrid(T=end - start, dt=dt, nt=nt, start=start)
    # validate time grid
    def _validate_time_grid(self):
        if self.time_grid is None:
            return
        start = self.time_grid.start
        end = self.time_grid.start + self.time_grid.T
        dt = self.time_grid.dt
        if start is None:
            log.error("time_grid.start is required")
        if end is None:
            log.error("time_grid.end is required")
        if dt is None:
            log.error("time_grid.dt is required")
        if start is not None and end is not None and end <= start:
            log.error("time_grid.end must be greater than time_grid.start")
        if dt is not None and dt <= 0:
            log.error("time_grid.dt must be > 0")
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
        if self.local_env_data is None:
            return
        if self.environment_config.is_explicit:
            self._validate_explicit_environment_data()
            return
        # validate common data
        self._validate_common_environment_data()
        # validate separate environments
        if self.env_model == "volcanic_rock":
            self._validate_volcanic_rock_environment_data()
        elif self.env_model == "surface_pond":
            self._validate_surface_pond_environment_data()
    # validate common local environment data
    def _validate_common_environment_data(self):
        if self.local_env_data.get("solvent_data") is not None:
            self._validate_liquid_level_params()
    # validate volcanic rock environment data
    def _validate_volcanic_rock_environment_data(self):
        missing = self._REQUIRED_ENV_VOLCROCK_KEYS - set(self.local_env_data.keys())
        if missing:
            log.error(f"Missing volcanic rock environment keys: {missing}")
        self._validate_temperature_pressure()
        n_pores = self.local_env_data.get("num_pores")
        if not isinstance(n_pores, int) or n_pores <= 0:
            log.error("number_pores must be a positive integer")
        # LOCAL VARIABLES
        pore_radius = self.local_env_data.get("pore_radius")
        pore_height = self.local_env_data.get("pore_height")
        min_distance = self.local_env_data.get("distance_neigh_pores")
        if pore_radius is None or pore_radius <= 0:
            log.error("pore_radius must be > 0")
        if pore_height is None or pore_height <= 0:
            log.error("pore_height must be > 0")
        if min_distance is None or min_distance <= 0:
            log.error("min_distance must be > 0")
    # validate surface pond environment data
    def _validate_surface_pond_environment_data(self):
        self._validate_temperature_pressure()
        pool_geometry = self.local_env_data.get("pool_geometry")
        if pool_geometry is None:
            return
        height = pool_geometry.get("height")
        surface_area_z0 = pool_geometry.get("surface_area_z0")
        if height is None or height <= 0:
            log.error("pool_geometry.height must be > 0")
        if surface_area_z0 is None or surface_area_z0 <= 0:
            log.error("pool_geometry.surface_area_z0 must be > 0")
    # validate single temperature / pressure environments
    def _validate_temperature_pressure(self):
        T = self.local_env_data.get("temperature")
        P = self.local_env_data.get("pressure")
        if T is None and not self.environment_config.uses_planetary_solver:
            log.error("Missing environment temperature")
        if P is None and not self.environment_config.uses_planetary_solver:
            log.error("Missing environment pressure")
        if T is not None and T <= 0:
            log.error("temperature must be > 0")
        if P is not None and P <= 0:
            log.error("pressure must be > 0")
    # validate explicit environment data
    def _validate_explicit_environment_data(self):
        if not self.local_env_data:
            log.error("Explicit environment source requires environment_data")
        if self.env_model == "hydro_vent":
            required_keys = {"T_hot", "T_cold", "pressure", "pH_hot", "pH_cold"}
        elif self.env_model == "surface_pond":
            required_keys = {"temperature", "pressure"}
        else:
            required_keys = set()
        missing = required_keys - set(self.local_env_data.keys())
        if missing:
            log.error(f"Missing explicit environment keys: {missing}")
    # validate liquid level parameters
    def _validate_liquid_level_params(self):
        params = self.local_env_data.get("liquid_level_params")
        if params is None:
            log.error("Missing liquid_level_params in environment data")
            return
        if params.base_level is None and not self.environment_config.uses_planetary_solver:
            log.error("liquid_level_params: 'base_level' is required")
