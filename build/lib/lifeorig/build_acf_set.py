import numpy as np
import random
from lifeorig.reaction_network import reaction_net_class
from lifeorig.read_input import p
#
# ACF set builder :
# make a sequence of catalysts set
# size = size_ACFS
def build_ACFS(size_ACFS, size_bpol, size_F, size_C):
    ACF_set = []
    i = 1
    while len(ACF_set) < size_ACFS:
        ACFS = reaction_net_class(size_bpol, size_F, size_C)
        # set up catalyst set
        size_X = ACFS.size_X
        catalyst_set = build_catalysts_set(size_X, size_C)
        ACFS.set_binary_polymer_model(catalyst_set)
        # prepare network plot
        file_name = p.working_dir+'/acs-' + str(i) + '.html'
        ACFS.show_network(file_name)
        # produce network genome
        ACFS.set_network_genome()
        ACF_set.append(ACFS)
        i += 1
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
def compute_hamming_distance(ACFS1, ACFS2):
    pass