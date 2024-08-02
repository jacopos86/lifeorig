from lifeorig.logging_module import log
from lifeorig.read_input import p
import numpy as np
import logging
from math import sin, exp
#
#  This module contains
#  routines for mutation rates
#

#
# compute Hamming distance 
# matrix
def compute_hamm_dist_matrix(ACF_set):
    N = len(ACF_set)
    HD_ij = np.zeros((N,N))
    # run over different CRS
    for i in range(N):
        g_1 = ACF_set[i].genome
        for j in range(N):
            g_2 = ACF_set[j].genome
            HD_ij[i,j] = compute_hamming_distance(g_1, g_2)
    return HD_ij

#
# compute Hamming distance
def compute_hamming_distance(ACFS1, ACFS2):
    N1 = len(ACFS1)
    N2 = len(ACFS2)
    if N1 != N2:
        log.error("\t len(S1) != len(S2)")
    N = N1
    d = 0
    for i in range(N):
        c1 = ACFS1[i]
        c2 = ACFS2[i]
        if c1 != c2:
            d += 1
    return d

#
# mutation site function

#
#   mutation rate class -> Q = I

class zero_mutation():
    def __init__(self, size):
        self.size = size
        self.Q = np.zeros((size, size))
        self.Q_oft = None
    def set_mut_matrix(self):
        self.Q = np.identity(self.size)
    # set time random matrix
    def set_constant_mutation_over_time(self, nt):
        self.Q_oft = np.zeros((self.size, self.size, nt))
        # iterate over t
        for t in range(nt):
            self.Q_oft[:,:,t] = self.Q[:,:]
    # show Q matrix
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
            
#
#    dist. mutation class
#

class dist_mutation():
    def __init__(self, size):
        self.size = size
        self.Q = np.zeros((size, size))
        self.Q_oft = None
    #
    # single site mutation
    def set_site_mut_oft(self, r, A, w, tau, t):
        # mu(t) = r + A sin(wt) * exp(-t/tau)
        # r is background mutation rate
        #
        nt = len(t)
        mu_oft = np.zeros(nt)
        # compute mu_oft
        for i in range(nt):
            mu_oft[i] = r + A * sin(w*t[i]) * exp(-t[i]/tau)
            if mu_oft[i] > 1.:
                mu_oft[i] = 1.
            elif mu_oft[i] < 0.:
                mu_oft[i] = 0.
        return mu_oft
    #
    # write mu_oft to file
    def write_muoft_to_file(self, t, mu_oft, out_file):
        nt = len(t)
        # open file
        f = open(out_file, 'w')
        for i in range(nt):
            f.write("%.10f          " % t[i] + "%.10f\n" % mu_oft[i])
        f.close()
    #
    # mutation matrix
    def set_mut_matrix(self, mu, n, d_H):
        # use the Eigen
        # definition
        # m_ij = mu^d_ij (1-mu)^(n-d_ij)
        # n = genome length
        # d_H : Hamming distance
        for i in range(self.size):
            for j in range(self.size):
                if i != j:
                    self.Q[i,j] = mu ** d_H[i,j] * (1. - mu) ** (n - d_H[i,j])
            self.Q[i,i] = 1. - sum(self.Q[i,:])
    #
    # mutation matrix of t
    def set_mut_matrix_oft(self, mu_oft, n, d_H, nt):
        self.Q_oft = np.zeros((self.size, self.size, nt))
        # use the Eigen
        # definition
        # m_ij = mu^d_ij (1-mu)^(n-d_ij)
        # n = genome length
        # d_H : Hamming distance
        for i in range(self.size):
            for j in range(self.size):
                if i != j:
                    for t in range(nt):
                        self.Q_oft[i,j,t] = mu_oft[t] ** d_H[i,j] * (1. - mu_oft[t]) ** (n - d_H[i,j])
                        #print(mu_oft[t] ** d_H[i,j], (1. - mu_oft[t]) ** (n - d_H[i,j]))
            for t in range(nt):
                self.Q_oft[i,i,t] = 1. - sum(self.Q_oft[i,:,t])
    #
    # show Q matrix
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