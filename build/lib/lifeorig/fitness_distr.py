import numpy as np
from lifeorig.read_input import p
from lifeorig.logging_module import log
# This class defines the fitness
# function
class fitness_distr():
    def __init__(self, n):
        # set fitness array
        self.fitness = np.zeros(n)
        self.fitness_oft = None
        # shape
        self.size = n
    def set_random_initial_fitness(self, s, M):
        np.random.seed(s)
        # set [0, M) random distribution
        self.fitness = np.random.rand(self.size)
        self.fitness[:] = self.fitness[:] * M
    def set_constant_fitness_over_time(self, nt):
        self.fitness_oft = np.zeros((self.size,nt))
        for t in range(nt):
            self.fitness_oft[:,t] = self.fitness[:]
    def show_fitness_distr(self):
        log.info("\t " + p.sep)
        log.info("\n")
        log.info("\t fitness distr : ")
        log.info("\n")
        line = ""
        for j in range(self.size):
            line += " {0:.3f}".format(self.fitness[j])
        log.info("\t " + line)
        log.info("\n")
        log.info("\t " + p.sep)