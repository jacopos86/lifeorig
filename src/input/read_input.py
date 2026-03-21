import json
import os
import numpy as np
from src.utilities.logging_module import log

#
#  parameters class

class parameters_class:
    ''' parameters class '''
    _ALLOWED_METABOLITE_TYPES = {"binary", "multi"}
    _ALLOWED_METABOLITE_DISTR = {"uniform", "length_decay"}
    def __init__(self):
        # work dir
        self.working_dir = None
        # n. protocell (initial)
        self.QSP_size = None
        # catalyst set parameters
        self.catalyst_set_params = None
        # rates distribution parameters
        self.rates_params = None
        # molecules data parameters
        self.metabolites_params = None
    def read_input_json(self, json_file):
        try:
            f = open(json_file)
        except:
            msg = "file: " + json_file + " not found"
            raise Exception(msg)
        data = json.load(f)
        f.close()
        # read input parameters
        if "working_dir" in data:
            self.working_dir = data["working_dir"]
            isExist = os.path.exists(self.working_dir)
            if not isExist:
                os.mkdir(self.working_dir)
        # num. individuals in QSP to average
        if "QSP_size" in data:
            self.QSP_size = data["QSP_size"]
        # metabolites data
        if "metabolites_data" in data:
            self.metabolites_params = data["metabolites_data"]
        # catalysts set size
        if "catalyst_set" in data:
            self.catalyst_set_params = data["catalyst_set"]
        #
        # mutation parameters
        #
        if "distribution_rates" in data:
            self.rates_params = data["distribution_rates"]
            print(self.rates_params)
        # time variables
        # size sample space
        if "evol_params" in data:
            self.evol_size = data["evol_params"]
    def validate(self):
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
            log.error.append(
                f"Invalid distr. type '{self.metabolites_params.get('metabolites_distr_type')}'. "
                f"Valid options: {sorted(self._ALLOWED_METABOLITE_DISTR)}"
            )

p = parameters_class()
p.sep = "*"*94