#  
#   main program :
#   life origin simulator
#
import numpy as np
import logging
from src.input.parser import parser
from src.input.read_input import p
from src.common.set_rndm_matrix import random_matrix
from src.mutations.mutation_rate import zero_mutation, dist_mutation
from src.QSP.fitness_distr import fitness_distr, fitness_distr_game_dyn
from src.QSP.quasi_species_solver import BuildQuasiSpeciesSolver
from src.chem_network.build_network import build_chem_networks
from src.chem_network.build_acfs_popul import build_ACFS_networks
from src.cell.build_QSP_list import set_up_empty_QSP_list
from src.molecules_dyn.gillespie_algo import chemical_kinetics_solver
from src.mutations.mutation_rate import compute_hamm_dist_matrix
from src.utilities.logging_module import log
from src.molecules.define_molecule_set import build_molecule_set

args = parser.parse_args()
calc_type = args.ct[0]
p.read_input_json(args.json_input[0])
p.validate()

log.info("\t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
log.info("\t ++++++                                                                                  ++++++")
log.info("\t ++++++                           LIFEORIG   CODE                                        ++++++")
log.info("\t ++++++                                                                                  ++++++")
log.info("\t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
log.info("\n")
log.info("\t " + p.sep)
log.info("\t CALCULATION TYPE : " + calc_type)
log.info("\t " + p.sep)
log.info("\n")

# test gillespie algo
if log.level <= logging.DEBUG:
    kin_solver = chemical_kinetics_solver()
    kin_solver.test()

# build chemical networks section

if calc_type == "set_initial_state":
    
    # set list molecular types

    X_set, X_set_map = build_molecule_set(p.metabolites_params)
    
    # set list of protocells

    n_protocells = p.QSP_size

    protocell_list = set_up_empty_QSP_list(n_protocells)

    protocell_list.set_initial_molecule_distr(X_set, p.metabolites_params)

    # set random networks -> one for each protocell

    build_chem_networks(protocell_list, X_set, p.bpol_strng_size, p.catalyst_set_params, p.rates_params)
    
# evolutionary section

elif calc_type == "evol":

    if p.fitness_eval != "compute":
        log.error("compute everything in evol")
    
    # first build sample space
    # from 1 to 10 of distance from
    # initial network

    for ic in range(1, p.n_acf_distr+1):
        
        ACF_distr = set_sample_space(p.evol_size, p.bpol_strng_size, p.size_F, p.size_C, ic)

        log.info("\t DISTRIBUTION LENGTH : " + str(len(ACF_distr)))

        # check if evolutionary game dynamics
        # has to be performed

        if p.EvolutionaryGameDyn:
            
            log.info("\t " + p.sep)
            log.info("\t PERFORM EVOLUTIONARY GAME DYNAMICS")
            log.info("\n")

            fitness_func = fitness_distr_game_dyn(len(ACF_distr))
            out_file = p.working_dir + "/" + str(ic) + "/a_ij.txt"
            fitness_func.set_fitness_distr(ACF_distr, out_file)

            log.info("\n")
            log.info("\t " + p.sep)
            log.info("\t PAYOFF MATRIX CALCULATION COMPLETED")
            log.info("\t " + p.sep)
            log.info("\n")

        else:
            
            log.info("\t " + p.sep)
            log.info("\t EXTRACT FITNESS FUNCTION")
            log.info("\n")

            # set fitness function

            fitness_func = fitness_distr(len(ACF_distr))
            fitness_func.set_fitness_distr(ACF_distr)
            fitness_func.show_fitness_distr()

            log.info("\n")
            log.info("\t " + p.sep)
            log.info("\t FITNESS CALCULATION COMPLETED")
            log.info("\t " + p.sep)
            log.info("\n")

        # initial distrib.

        x0 = np.zeros(len(ACF_distr))
        x0[0] = 1.
        x0[:] = x0[:] / np.sum(x0)
        assert np.abs(sum(x0)-1.) < 1.E-7
        log.info("\t x0 = " + str(x0[:min(50,len(ACF_distr))]))
        log.info("\n")
        log.info("\t " + p.sep)

        # time variables

        log.info("\n")
        log.info("\t dt : " + str(p.dt))
        log.info("\t T  : " + str(p.T))
        log.info("\t " + p.sep)

        # set mutation matrix

        if p.mutation_typ == "random":
            mutation_obj = random_matrix(len(ACF_distr))
            mutation_obj.set_rand_matrix(p.seed, p.r_mut, p.sig_mut)
            mutation_obj.normalize_matrix()
            mutation_obj.show_Q_matrix()
        elif p.mutation_typ == "zero":
            mutation_obj = zero_mutation(len(ACF_distr))
            mutation_obj.set_mut_matrix()
            mutation_obj.show_Q_matrix()
        elif p.mutation_typ == "dist":
            mutation_obj = dist_mutation(len(ACF_distr))
            n = len(ACF_distr[0].genome)
            # compute Hamming distance
            HDij = compute_hamm_dist_matrix(ACF_distr)
            n = min(n, HDij.max())
            mutation_obj.set_mut_matrix(p.r_mut, n, HDij)
            mutation_obj.show_Q_matrix()

        # set quasi species solver

        solver = BuildQuasiSpeciesSolver(len(ACF_distr), p.dt, p.T)
        xt = solver.solve(x0, fitness_func, mutation_obj, ACF_distr)
        out_file = p.working_dir + "/" + str(ic) + "/x_oft.txt"
        out_file2= p.working_dir + "/" + str(ic) + "/Q_oft.txt"
        solver.save_data(xt, out_file, out_file2)
    
log.info("\t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
log.info("\t ++++++                                                                                  ++++++")
log.info("\t ++++++                 LIFEORIG   CODE    EXECUTION   COMPLETE                          ++++++")
log.info("\t ++++++                                                                                  ++++++")
log.info("\t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
log.info("\n")