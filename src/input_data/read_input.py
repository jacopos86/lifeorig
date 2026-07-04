import os
from src.utilities.logging_module import log
from src.input_data.core import AbstractInput
from src.input_data.input_configs import (
    ChemicalNetworkInputConfig,
    EnvironmentInputConfig,
    EvolutionInputConfig,
    MoleculeInputConfig,
)
from src.input_data.environment_input import EnvironmentInputBuilder
from src.input_data.planet_input import PlanetInputBuilder
from src.environment.set_planet_environment import derive_planet_env_data
from src.chemical_env.chemical_environment import set_ChemEnvInput, SimulateChemEnv

#
#  parameters class

class parameters_class(AbstractInput, PlanetInputBuilder, EnvironmentInputBuilder):
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
        self.planetary_data = {}
        self.chem_env_params = None
        # stellar data
        self.stellar_data = None
        # evolution data
        self.evolution_config = None
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
    # set local environment data and derived planetary forcing
    def set_local_environment(self):
        self.env_model = self.environment_config.env_type
        self.local_env_data = self._build_local_env_data()
        if self.local_env_data is None:
            return
        if self.environment_config.is_explicit:
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
    #
    #    validation section
    #
    def _validate(self):
        # optional: check metabolites_parameters
        required_keys = ["type", "initial_population_molecules"]
        if self.metabolites_params.get("type") in {"binary", "multi"}:
            required_keys.extend(["pol_strng_maxsize", "metabolites_distr_type"])
        if self.metabolites_params.get("type") == "reference_file":
            required_keys.append("reaction_file")
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
        if self.local_env_data is None:
            return
        if self.environment_config.is_explicit:
            self._validate_explicit_environment_data()
            return
        if self.env_model == "volcanic_rock":
            # check keys
            missing = self._REQUIRED_ENV_VOLCROCK_KEYS - set(self.local_env_data.keys())
            if missing:
                log.error(f"Missing environment keys: {missing}")
            # number of pores
            n_pores = self.local_env_data["num_pores"]
            if not isinstance(n_pores, int) or n_pores <= 0:
                log.error("number_pores must be a positive integer")
            # pore geometry
            pore_radius = self.local_env_data["pore_radius"]
            pore_height = self.local_env_data["pore_height"]
            if pore_radius <= 0:
                log.error("pore_radius must be > 0")
            if pore_height <= 0:
                log.error("pore_height must be > 0")
            # distance between closest neighbors
            min_distance = self.local_env_data["distance_neigh_pores"]
            if min_distance <= 0:
                log.error("min_distance must be > 0")
        # ----------------------------
        # temperature / pressure
        # ----------------------------
        T = self.local_env_data["temperature"]
        P = self.local_env_data["pressure"]
        if T <= 0:
            log.error("temperature must be > 0")
        if P <= 0:
            log.error("pressure must be > 0")
        self._validate_liquid_level_params()
    # validate explicit environment data
    def _validate_explicit_environment_data(self):
        if not self.local_env_data:
            log.error("Explicit environment source requires environment_data")
        if self.env_model == "hydro_vent":
            required_keys = {"T_hot", "T_cold", "pressure", "pH_hot", "pH_cold"}
        elif self.env_model == "surface_pond":
            required_keys = {"temperature", "pressure", "pH"}
        else:
            required_keys = set()
        missing = required_keys - set(self.local_env_data.keys())
        if missing:
            log.error(f"Missing explicit environment keys: {missing}")
    # validate liquid level parameters
    def _validate_liquid_level_params(self):
        solvent_data = self.local_env_data.get("solvent_data")
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
