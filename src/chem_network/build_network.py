import logging
from src.chem_network.reaction_network import ReactionNetwork
from src.input.read_input import p
from src.utilities.logging_module import log
from src.catalysts.catalysts_set import set_catalysts_prob_distr
#
# ACF set builder :
# make a sequence of catalysts set
# size = size_ACFS
def build_chem_networks(protocell_list, X_set, size_bpol, catalyst_set_params, rates_params):
    # iterate over protocells
    for protocell in protocell_list:
        print(protocell)
        Rnet = ReactionNetwork(size_bpol)
        # set up catalyst set
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