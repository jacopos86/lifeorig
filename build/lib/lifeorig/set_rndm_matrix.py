import numpy as np
from lifeorig.logging_module import log
from lifeorig.read_input import p
import logging
# This module set up a random matrix
# conditions to be satisfied :
# 1) p_ij >= 0
# 2) sum_j p_ij = 1
# 3) at least one element in each column .ne. 0
class random_matrix():
    def __init__(self, size):
        self.size = size
        self.Q = np.zeros((size,size))
        self.Q_oft = None
    def set_rand_matrix(self, s):
        np.random.seed(s)
        # uniform distrib. [0,1)
        self.Q = np.random.rand(self.size,self.size)
    def normalize_matrix(self):
        for i in range(self.size):
            s = sum(self.Q[i,:])
            self.Q[i,:] = self.Q[i,:] / s
    # set random matrix with pointwise
    # mutations
    def set_rand_matrix_pwmut(self):
        pass
    # set time random matrix
    def set_constant_mutation_over_time(self, nt):
        self.Q_oft = np.zeros((self.size, self.size, nt))
        # iterate over t
        for t in range(nt):
            self.Q_oft[:,:,t] = self.Q[:,:]
    def show_Q_matrix(self):
        log.info("\t " + p.sep)
        log.info("\n")
        log.info("\t Q matrix : ")
        log.info("\n")
        for i in range(self.size):
            line = ""
            for j in range(self.size):
                line += "  {0:.3f}".format(self.Q[i,j])
            log.info("\t " + line)
        log.info("\n")
        log.info("\t " + p.sep)
        # debug section
        if log.level <= logging.DEBUG:
            log.info("\n")
            for i in range(self.size):
                log.debug("\t row {}".format(i+1) + " - \t sum : {0:.3f}".format(sum(self.Q[i,:])))
            log.info("\n")
            log.info("\t " + p.sep)