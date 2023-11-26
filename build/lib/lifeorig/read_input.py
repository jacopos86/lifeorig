import json
import os
#
#  parameters class

class parameters_class:

    def __init__(self):
        self.GA_type = None

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
        if "GA_type" in data:
            self.GA_type = data["GA_type"]
        if "num_pop" in data:
            self.n_pop = data["num_pop"]
        if "n_bits" in data:
            self.n_bits = data["n_bits"]
        if "r_mut" in data:
            self.r_mut = data["r_mut"]
        if "r_cross" in data:
            self.r_cross = data["r_cross"]
        if "max_iter" in data:
            self.n_max_iter = data["max_iter"]

p = parameters_class()