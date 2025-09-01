import logging
from reaction_network import reaction_net_class
from read_input import p
from logging_module import log
from catalysts_set import build_catalysts_rr_set, set_catalysts_prob_distr
#
# ACF set builder :
# make a sequence of catalysts set
# size = size_ACFS
def build_chem_networks(num_ACFS, size_bpol, size_F, size_C, size_C_intern, ratio_C_ACFset, nQSP_indiv):
    # open data file
    file_name = p.working_dir + "/ACF_data.txt"
    f = open(file_name, 'a')
    ACF_set_list = []
    i = 1
    while len(ACF_set_list) < num_ACFS:
        ACFS = reaction_net_class(size_bpol, size_F, size_C)
        # set up catalyst set
        size_X = ACFS.size_X
        print(size_X, size_C, ratio_C_ACFset[i-1])
        catalyst_rr_set = build_catalysts_rr_set(size_C, size_C_intern[i-1])
        catalyst_prob_distr = set_catalysts_prob_distr(size_C, size_C_intern[i-1], ratio_C_ACFset[i-1])
        exit()
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