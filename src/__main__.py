#  
#   main program :
#   life origin simulator
#
import numpy as np
import logging
from src.input_data.parser import parser
from src.mutations.mutation_rate import zero_mutation, dist_mutation
from src.QSP.fitness_distr import fitness_distr, fitness_distr_game_dyn
from src.QSP.quasi_species_solver import BuildQuasiSpeciesSolver
from src.chem_network.build_network import build_chem_networks
from src.chem_network.build_acfs_popul import build_ACFS_networks
from src.cell.build_QSP_list import set_up_empty_QSP_list
from src.molecules_dyn.gillespie_algo import chemical_kinetics_solver
from src.mutations.mutation_rate import compute_hamm_dist_matrix
from src.utilities.logging_module import log
from src.metabolites.metabolite_builder import build_metabolites
from src.catalysts.catalysts_set import build_catalyst_set
from src.environment.setup_environment import set_simulation_environment
from src.input_data.read_input import parameters_class
from src.environment.setup_environment import setup_environment
from src.environment.planetary_driver import planetary_solver_driver
from src.network_generation.reaction_database_driver import reaction_database_driver
from src.network_generation.reaction_mysql_db import open_reaction_database, log_species_summary
from src.chemical_types.define_molecule_set import build_molecular_species_set

args = parser.parse_args()
calc_type = args.ct[0]
p = parameters_class()
p.read_input_json(args.json_input[0])

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

    if p.environment_config.uses_planetary_solver:
        log.info("\t " + p.sep)
        log.info("\t STARTING PLANETARY SOLVER")
        planet_params = planetary_solver_driver(p)
        p.planetary_data = planet_params
        log.info("\n")
    else:
        log.info("\t SKIPPING PLANETARY SOLVER: environment_source is explicit")

    # reaction database driver

    reaction_source_files = reaction_database_driver(p)
    db = open_reaction_database()
    try:
        log_species_summary(db, reaction_source_files)
    finally:
        db.close()
    
    # set up full chemical set

    species_set = build_molecular_species_set(p, reaction_source_files=reaction_source_files)
    for x in species_set.molecules:
        print(x)
    print("\n")
    for x in species_set.templates:
        print(x)
    print("\n")
    print(species_set.minerals)
    exit()
    # set up local environment

    log.info("\t " + p.sep)
    log.info("\t SETTING UP ENVIRONMENT")
    env_setup = setup_environment(p, species_set)
    log.info("\n")
    exit()
    
    # set list molecular types
    
    X_set, X_set_map, X_init = build_metabolites(p.metabolites_params, species_set)
    exit()
    # build catalysts set: Y set

    Y_set = build_catalyst_set(X_set, p.catalyst_set_params)
    
    # set list of protocells

    n_protocells = p.QSP_size

    protocell_list = set_up_empty_QSP_list(n_protocells, p.protocell_info)

    # create simulation environment
    
    simul_env = set_simulation_environment(
        protocell_list,
        env_type=p.env_model,
        env_data=p.env_data
    )

    # set random networks -> one for each protocell

    build_chem_networks(simul_env, X_set, Y_set, p.metabolites_params, p.rates_params)
    
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
