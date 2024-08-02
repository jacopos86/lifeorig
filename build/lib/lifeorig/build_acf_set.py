import numpy as np
import random
import logging
from lifeorig.reaction_network import reaction_net_class
from lifeorig.read_input import p
from lifeorig.logging_module import log
#
# ACF set builder :
# make a sequence of catalysts set
# size = size_ACFS
def build_ACFS(size_ACFS, size_bpol, size_F, size_C):
    # open data file
    file_name = p.working_dir + "/ACF_data.txt"
    f = open(file_name, 'a')
    ACF_set = []
    i = 1
    while len(ACF_set) < size_ACFS:
        ACFS = reaction_net_class(size_bpol, size_F, size_C)
        # set up catalyst set
        size_X = ACFS.size_X
        catalyst_set = build_catalysts_set(size_X, size_C)
        ACFS.set_binary_polymer_model(catalyst_set)
        ACFS.find_ACF_subset()
        # prepare network plot
        if log.level == logging.DEBUG:
            file_name = p.working_dir+'/acs-' + str(i) + '.html'
            ACFS.show_network(file_name)
        # produce network genome
        ACFS.set_network_genome()
        # here we solve the kinetic model
        # multiple times -> average different final
        # configurations
        if p.fitness_eval == "compute":
            ACFS.set_chemical_kinetics_solver()
            # append to file
            f.write("%s          " % ACFS.genome + "%d          " % len(ACFS.ACF_set) + "%.7f\n" % ACFS.fitness)
            # append data to file
        ACF_set.append(ACFS)
        log.info("\n")
        log.info("\t " + p.sep)
        log.warning("\t ACF SET : " + str(i))
        log.info("\t " + p.sep)
        log.info("\n")
        i += 1
    f.close()
    return ACF_set
#
# build catalyst set
# 
def build_catalysts_set(size_X, size_C):
    catalyst_set = []
    while len(catalyst_set) < size_C:
        i = random.choice(np.arange(1, size_X+1, 1))
        if i not in catalyst_set:
            catalyst_set.append(i)
    return catalyst_set