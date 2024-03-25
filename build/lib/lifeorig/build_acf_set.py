import numpy as np
import random
import logging
from lifeorig.reaction_network import reaction_net_class
from lifeorig.read_input import p
from lifeorig.logging_module import log
#
# ACF set builder :
# make a sequence of catalysts set
# size = size_ACFS
def build_ACFS(size_ACFS, size_bpol, size_F, size_C):
    # open data file
    file_name = p.working_dir + "/ACF_data.txt"
    f = open(file_name, 'a')
    ACF_set = []
    i = 1
    while len(ACF_set) < size_ACFS:
        ACFS = reaction_net_class(size_bpol, size_F, size_C)
        # set up catalyst set
        size_X = ACFS.size_X
        catalyst_set = build_catalysts_set(size_X, size_C)
        ACFS.set_binary_polymer_model(catalyst_set)
        ACFS.find_ACF_subset()
        # prepare network plot
        if log.level == logging.DEBUG:
            file_name = p.working_dir+'/acs-' + str(i) + '.html'
            ACFS.show_network(file_name)
        # produce network genome
        ACFS.set_network_genome()
        # here we solve the kinetic model
        # multiple times -> average different final
        # configurations
        ACFS.set_chemical_kinetics_solver()
        ACF_set.append(ACFS)
        # append to file
        f.write("%s          " % ACFS.genome + "%d          " % len(ACFS.ACF_set) + "%.7f\n" % ACFS.fitness)
        # append data to file
        i += 1
    f.close()
    return ACF_set
#
# build catalyst set
# 
def build_catalysts_set(size_X, size_C):
    catalyst_set = []
    while len(catalyst_set) < size_C:
        i = random.choice(np.arange(1, size_X+1, 1))
        if i not in catalyst_set:
            catalyst_set.append(i)
    return catalyst_set
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