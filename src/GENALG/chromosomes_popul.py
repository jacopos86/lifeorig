from lifeorig.read_input import p
import numpy as np
from numpy.random import randint
from numpy.random import rand
from lifeorig.objective_func import objective
from lifeorig.logging_module import log
#
#  interface classical / quantum
#  GA
class chromosomes_pop_class:
    def __init__(self):
        # n. chromosomes in population
        self.n_bits = p.n_bits
        self.n_pop = p.n_pop
    # set up population
    def set_initial_population(self):
        self.pop = [randint(0, 2, self.n_bits).tolist() for _ in range(self.n_pop)]
    def update_population(self, new_pop):
        self.pop = []
        self.pop = new_pop
    def mutation(self, bit_string, r_mut):
        for i in range(len(bit_string)):
            # check for mutation
            if rand() < r_mut:
                # flip the bit
                bit_string[i] = 1 - bit_string[i]
        return bit_string
    def crossover(self, p1, p2, r_cross):
        # by default children are simple copies
        c1, c2 = p1.copy(), p2.copy()
        # check recombination
        if rand() < r_cross:
            # select crossover pt
            pt = randint(1, len(p1)-2)
            # crossover
            c1 = p1[:pt] + p2[pt:]
            c2 = p2[:pt] + p1[pt:]
        return [c1, c2]
    # produce next gen.
    def produce_next_gen(self, selected_pop):
        children = []
        for i in range(0, self.n_pop, 2):
            # get the parents
            p1, p2 = selected_pop[i], selected_pop[i+1]
            # crossover + mutation
            for c in self.crossover(p1, p2, p.r_cross):
                # mutation
                c = self.mutation(c, p.r_mut)
                children.append(c)
        self.update_population(children)
    # evaluate the score of all candidates
    # in the population
    def compute_scores(self):
        self.scores = [objective(c) for c in self.pop]
    def print_avg_score(self):
        avg_score = np.array(self.scores).mean()
        log.info("avg. score : " + str(avg_score))