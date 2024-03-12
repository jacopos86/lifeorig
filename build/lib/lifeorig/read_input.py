import json
import os
import numpy as np
#
#  parameters class

class parameters_class:

    def __init__(self):
        # max. possible fitness
        self.max_fitness = 1.
        # random seed
        self.seed = 1

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
        # genome size
        if "ACFS_size" in data:
            self.ACFS_size = data["ACFS_size"]
        # max fitness
        if "max_fitness" in data:
            self.max_fitness = data["max_fitness"]
        # food set size
        if "food_set_size" in data:
            self.size_F = data["food_set_size"]
        # catalysts set size
        if "catalyst_set_size" in data:
            self.size_C = data["catalyst_set_size"]
        # molecule set size
        if "bpol_strng_size" in data:
            self.bpol_strng_size = data["bpol_strng_size"]
        # total molecules population
        if "total_population_molecules" in data:
            self.n0 = data["total_population_molecules"]
        # time variables
        if "T" in data:
            self.T = data["T"]
        if "dt" in data:
            self.dt = data["dt"] 
        if "r_mut" in data:
            self.r_mut = data["r_mut"]
        if "r_cross" in data:
            self.r_cross = data["r_cross"]
        if "max_iter" in data:
            self.n_max_iter = data["max_iter"]
        if "random_seed" in data:
            self.seed = data["random_seed"]

p = parameters_class()
p.sep = "*"*94