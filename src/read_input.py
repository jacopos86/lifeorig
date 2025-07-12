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
        # n. config
        self.nconfig = 1
        # n. ACF distr
        self.n_acf_distr = 1
        # set evolutionary game dyn.
        self.EvolutionaryGameDyn = False

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
        # evolutionary game dyn.
        if "EvolutionaryGameDyn" in data:
            self.EvolutionaryGameDyn = data["EvolutionaryGameDyn"]
        # size initial ACF set
        if "num_ACFS" in data:
            self.num_ACFS = data["num_ACFS"]
        # num. individuals in QSP to average
        if "n_QSP_indiv" in data:
            self.QSP_indiv = data["n_QSP_indiv"]
        # size sample space
        if "evol_size" in data:
            self.evol_size = data["evol_size"]
        # fitness eval. model
        if "fitness_eval" in data:
            self.fitness_eval = data["fitness_eval"]
        # max fitness
        if "max_fitness" in data:
            self.max_fitness = data["max_fitness"]
        # food set size
        if "food_set_size" in data:
            self.size_F = data["food_set_size"]
        # catalysts set size
        if "catalyst_set_size" in data:
            self.size_C = data["catalyst_set_size"]
        # internal network catalyst set size
        if "intern_catalyst_set_size" in data:
            self.intern_C_size = data["intern_catalyst_set_size"]
        if "ext_catalyst_concentr" in data:
            self.ext_catalyst_concentr = data["ext_catalyst_concentr"]
        # catalysts distribution
        if "ACFS_catalyst_prob_distr" in data:
            self.ACFS_catalyst_prob_distr = data["ACFS_catalyst_prob_distr"]
        if "ext_catalyst_prob_distr" in data:
            self.ext_catalyst_prob_distr = data["ext_catalyst_prob_distr"]
        if "ratio_C_ACF_set" in data:
            self.ratio_C_ACFset = data["ratio_C_ACF_set"]
        # molecule set size
        if "bpol_strng_size" in data:
            self.bpol_strng_size = data["bpol_strng_size"]
        # total molecules population
        if "total_population_molecules" in data:
            self.n0 = data["total_population_molecules"]
        # target molecules
        if "target_molecules" in data:
            self.target_molecules = data["target_molecules"]
        # n. config. to average 
        # kinetic simulation
        if "number_config" in data:
            self.nconfig = data["number_config"]
        # n. ACF distr.
        if "number_ACF" in data:
            self.n_acf_distr = data["number_ACF"]
        # mutation
        if "mutation" in data:
            self.mutation_typ = data["mutation"]
        # time variables
        if "T" in data:
            self.T = data["T"]
        if "dt" in data:
            self.dt = data["dt"]
        # mutation parameters
        #
        if "r_mut" in data:
            self.r_mut = data["r_mut"]
        if "sig_mut" in data:
            self.sig_mut = data["sig_mut"]
        if "A_mut" in data:
            self.A_mut = data["A_mut"]
        if "w_mut" in data:
            self.w_mut = data["w_mut"]
        if "tau_mut" in data:
            self.tau_mut = data["tau_mut"]
        # random mutation seed
        if "random_seed" in data:
            self.seed = data["random_seed"]
            
    def read_NN_params(self, json_file):
        # read file
        try:
            f = open(json_file)
        except:
            msg = "file: " + json_file + " not found"
            raise Exception(msg)
        data = json.load(f)
        f.close()
        # read input parameters
        if 'NN_model' in data:
            self.NN_model = data["NN_model"]
        if 'NN_parameters' in data:
            self.NN_parameters = data["NN_parameters"]

p = parameters_class()
p.sep = "*"*94