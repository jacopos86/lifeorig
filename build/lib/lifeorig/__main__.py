#  
#   main program :
#   life origin simulator
#
from lifeorig.parser import parser
from lifeorig.read_input import p
from lifeorig.set_rndm_matrix import random_matrix
from lifeorig.fitness_distr import fitness_distr
from lifeorig.quasi_species_solver import QuasiSpeciesSolver
from lifeorig.build_acf_set import build_ACFS
from lifeorig.logging_module import log

args = parser.parse_args()
p.read_input_json(args.json_input)

# set up GA

log.info("\t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
log.info("\t ++++++                                                                                  ++++++")
log.info("\t ++++++                           LIFEORIG   CODE                                        ++++++")
log.info("\t ++++++                                                                                  ++++++")
log.info("\t ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
log.info("\n")

#chromosomes = chromosomes_pop_class()
#chromosomes.set_initial_population()
#chromosomes.compute_scores()
#chromosomes.print_avg_score()

#for iter in range(p.n_max_iter):
#    selected = [tournament_selection(chromosomes) for _ in range(chromosomes.n_pop)]
#    chromosomes.produce_next_gen(selected)
#    chromosomes.compute_scores()
#    chromosomes.print_avg_score()

# set random matrix

size = p.ACFS_size
mutation_obj = random_matrix(size)
mutation_obj.set_rand_matrix(p.seed)
mutation_obj.normalize_matrix()
mutation_obj.show_Q_matrix()

# set fitness function

fitness_func = fitness_distr(size)
fitness_func.set_random_initial_fitness(p.seed, p.max_fitness)
fitness_func.show_fitness_distr()

# initial distrib.

log.info("\t x0 = " + str(p.x0))
log.info("\n")
log.info("\t " + p.sep)

# time variables

log.info("\n")
log.info("\t dt : " + str(p.dt))
log.info("\t T  : " + str(p.T))
log.info("\t " + p.sep)

# set quasi species solver

solver = QuasiSpeciesSolver(size, p.dt, p.T)
solver.solve(p.x0, fitness_func, mutation_obj)

# set random networks -> one for each species

ACF_set = build_ACFS(size, p.bpol_strng_size, p.size_F, p.size_C)