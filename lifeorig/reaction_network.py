import math
import random
import numpy as np
from lifeorig.logging_module import log
#  class describing the
#  reaction network -> we use a binary polymer model
class reaction_net_class:
    def __init__(self, size_F, size_C, size_X):
        # food set size
        self.size_F = size_F
        # catalysts size
        self.size_C = size_C
        self.catalyst_set = []
        # size of molecules set
        self.size_X = size_X
        # reactions set
        self.set_reactions = [{}]
    def set_binary_polymer_model(self):
        # max. size string
        # polymer model
        self.strng_size = math.log2(self.size_X)
        # n. food set bits
        self.n_F_bits = math.log2(self.size_F)
        log.info("\t max. string size : " + str(self.strng_size))
        # build the catalysts set (C)
        self.build_catalysts_set()
        # build set of reactions
        #self.build_reactions_set()
    # catalysts set
    def build_catalysts_set(self):
        while len(self.catalyst_set) <= self.size_C:
            i = random.choice(np.arange(0, self.size_X, 1))
            if i not in self.catalyst_set:
                self.catalyst_set.append(i)
    # reaction set building
    # method
    def build_reactions_set(self):
        # list of dictionaries with keys
        # for ligation
        # 1) 'r1' : first reactant
        # 2) 'r2' : second reactant
        # 3) 'c'  : catalyst
        # 4) 'p'  : product
        for r1 in range(self.size_X):
            x1 = bin(r1)
            for r2 in range(r1, self.size_X):
                x2 = bin(r2)
                print(x1, x2)
        # for cleavage
        # 1) 'r'  : reactant
        # 2) 'c'  : catalyst
        # 3) 'p1' : product 1
        # 4) 'p2' : product 2