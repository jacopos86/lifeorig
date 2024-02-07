import numpy as np
from lifeorig.logging_module import log
#  class describing the
#  reaction network
class reaction_net_class:
    def __init__(self, nn, nl):
        # n. nodes
        self.n_nodes = nn
        # n. links
        self.n_links = nl
    def build_rand_network(self):
        # build nodes list
        self.nodes = np.arange(1, self.n_nodes+1, 1)
        log.info("\t nodes list : " + str(self.nodes))