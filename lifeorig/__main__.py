#  
#   main program :
#   life origin simulator
#
import numpy as np
import random
from lifeorig.parser import parser
from lifeorig.read_input import p
from lifeorig.set_rndm_matrix import random_matrix
from lifeorig.fitness_distr import fitness_distr
from lifeorig.quasi_species_solver import QuasiSpeciesSolver
from lifeorig.build_acf_set import build_ACFS, compute_hamm_dist_matrix
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

x0 = np.zeros(p.ACFS_size)
for i in range(p.ACFS_size):
    x0[i] = random.uniform(0., 1.)
x0[:] = x0[:] / np.sum(x0)
assert np.abs(sum(x0)-1.) < 1.E-7
log.info("\t x0 = " + str(x0[:50]))
log.info("\n")
log.info("\t " + p.sep)

# time variables

log.info("\n")
log.info("\t dt : " + str(p.dt))
log.info("\t T  : " + str(p.T))
log.info("\t " + p.sep)

# set quasi species solver

solver = QuasiSpeciesSolver(size, p.dt, p.T)
solver.solve(x0, fitness_func, mutation_obj)

# set random networks -> one for each species

ACF_set = build_ACFS(size, p.bpol_strng_size, p.size_F, p.size_C)
HDij = compute_hamm_dist_matrix(ACF_set)
ACF_set[0].find_ACF_subset()