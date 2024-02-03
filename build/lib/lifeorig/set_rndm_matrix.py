import numpy as np
from lifeorig.logging_module import log
# This module set up a random matrix
# conditions to be satisfied :
# 1) p_ij >= 0
# 2) sum_j p_ij = 1
# 3) at least one element in each column .ne. 0
class random_matrix():
    def __init__(self, size):
        self.size = size
        self.Q = np.zeros((size,size))
    def set_rand_matrix(self):
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
    def show_Q_matrix(self):
        for i in range(self.size):
            for j in range(self.size):
                log.info(" " + str(self.Q[i,j]))
            log.info("\n")
            