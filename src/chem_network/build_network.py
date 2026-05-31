import logging
from src.input_data.read_input import p
from src.utilities.logging_module import log
from src.catalysts.catalysts_set import set_catalysts_prob_distr

#
# ACF set builder :
# make a sequence of chem. networks for each protocell
#

def build_chem_networks(simul_env, X_set, Y_set, metabolites_params, rates_params):
    '''
    Build set of chemical networks for each protocell
    '''
    # set initial metabolites distribution
    metabol_distr = simul_env.protocell_list.set_initial_metabolite_distr(X_set, metabolites_params)
    # iterate over protocells
    for protocell in simul_env.protocell_list:
        # --------------------------------------
        # 1)  set up the reaction network
        # --------------------------------------
        protocell.set_reaction_network(metabolites_params, X_set)
        # --------------------------------------
        # 2) define metabolites popul.
        # --------------------------------------
        protocell.set_initial_metabolite_population(metabolites_params, X_set, metabol_distr)
        protocell.set_initial_protocell_state()
        # --------------------------------------
        # 3)  set up catalytic state of set
        # --------------------------------------
        protocell.set_catalytic_state(X_set, Y_set)
        exit()
        catalyst_rr_set = build_catalysts_rr_set(size_C, size_C_intern[i-1])
        catalyst_prob_distr = set_catalysts_prob_distr(size_C, size_C_intern[i-1], ratio_C_ACFset[i-1])
        protocell.info()
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