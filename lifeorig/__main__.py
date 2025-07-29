#  
#   main program :
#   life origin simulator
#
import numpy as np
import logging
from lifeorig.parser import parser
from lifeorig.read_input import p
from lifeorig.set_rndm_matrix import random_matrix
from lifeorig.mutation_rate import zero_mutation, dist_mutation
from lifeorig.fitness_distr import fitness_distr, fitness_distr_game_dyn
from lifeorig.quasi_species_solver import BuildQuasiSpeciesSolver
from lifeorig.build_acf_set import build_ACFS
from lifeorig.build_sample_space import set_sample_space
from lifeorig.gillespie_algo import chemical_kinetics_solver
from lifeorig.mutation_rate import compute_hamm_dist_matrix
#from lifeorig.network_predict import generate_predictor, build_Xy_data
from lifeorig.logging_module import log

args = parser.parse_args()
calc_type = args.ct[0]
p.read_input_json(args.json_input[0])
if calc_type == "ml":
    p.read_NN_params(args.json_NN_input[0])

# set up GA

log.info("\t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
log.info("\t ++++++                                                                                  ++++++")
log.info("\t ++++++                           LIFEORIG   CODE                                        ++++++")
log.info("\t ++++++                                                                                  ++++++")
log.info("\t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
log.info("\n")
print(calc_type)

# test gillespie algo
if log.level <= logging.DEBUG:
    kin_solver = chemical_kinetics_solver()
    kin_solver.test()

if calc_type == "acf_set":
    
    size = p.QSP_size
    
    # set random networks -> one for each species

    ACF_set = build_ACFS(size, p.bpol_strng_size, p.size_F, p.size_C)

elif calc_type == "ml":
    
    size = p.ACFS_size
    
    # set random networks -> one for each species

    # here read ACF_set
    
    #   ACF_set = build_ACFS(size, p.bpol_strng_size, p.size_F, p.size_C)
    
    # read fitness distrib.
    
    # learn fitness model
    
    # genome_size = len(ACF_set[0].genome)
    # fitness_model = generate_predictor(p.NN_model, genome_size)
    # fitness_model.set_model(p.NN_parameters)
    # X, y = build_Xy_data(ACF_set)
    # X_test, y_test = fitness_model.fit(p.NN_parameters, X, y)
    # log.info("\t NN TEST SCORE : " + str(fitness_model.get_score(X_test, y_test)))
    
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