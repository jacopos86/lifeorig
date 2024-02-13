import math
import random
import numpy as np
import logging
from lifeorig.logging_module import log
from lifeorig.read_input import p
#  class describing the
#  reaction network -> we use a binary polymer model
class reaction_net_class:
    def __init__(self, size_bpol, size_F, size_C):
        # max. size string
        # polymer model
        self.strng_size = size_bpol
        # food set size
        self.size_F = size_F
        # catalysts size
        self.size_C = size_C
        self.catalyst_set = []
        # size of molecules set
        self.size_X = 2 ** size_bpol - 1
        # reactions set
        self.ligand_reactions = []
        self.cleavage_reactions = []
    def set_binary_polymer_model(self):
        # n. food set bits
        self.n_F_bits = math.log2(self.size_F)
        log.info("\t max. string size : " + str(self.strng_size))
        # build the catalysts set (C)
        self.build_catalysts_set()
        # build set of reactions
        self.build_reactions_set()
    # catalysts set
    def build_catalysts_set(self):
        while len(self.catalyst_set) <= self.size_C:
            i = random.choice(np.arange(1, self.size_X+1, 1))
            if i not in self.catalyst_set:
                self.catalyst_set.append(i)
        log.info("\n")
        log.info("\t C = " + str(self.catalyst_set))
        log.info("\t " + p.sep)
    # binary string -> decimal number
    def convert_binstr_to_dec(self, y):
        y = y.zfill(self.strng_size)
        n = self.strng_size-1
        r = 0
        for c in y:
            r += int(c) * 2 ** n
            n = n-1
        return int(r)
    # split string into
    # all possible substrings
    def split_strng_to_substrng(self, strng):
        substrngs = {strng[a:a+k] for k in range(1,1+len(strng)) 
            for a in range(1+len(strng)-k)}
        substr_list = []
        for sstr in substrngs:
            if len(sstr) != len(strng):
                substr_list.append(sstr)
        return substr_list
    # reaction set building
    # method
    def build_reactions_set(self):
        # list of dictionaries with keys
        # for ligation
        # 1) 'r1' : first reactant
        # 2) 'r2' : second reactant
        # 3) 'c'  : catalyst
        # 4) 'p'  : product
        for r1 in range(1, self.size_X+1):
            x1 = bin(r1)
            x1 = x1[2:]
            n1 = len(x1)
            for r2 in range(r1, self.size_X+1):
                x2 = bin(r2)
                x2 = x2[2:]
                n2 = len(x2)
                if n1 + n2 <= self.strng_size:
                    y = x1 + x2
                    # convert to decimal
                    pr= self.convert_binstr_to_dec(y)
                    # set dictionary
                    react = {}
                    react['r1'] = x1.rjust(self.strng_size, '.')
                    react['r2'] = x2.rjust(self.strng_size, '.')
                    react['p']  = y.rjust(self.strng_size, '.')
                    # select catalyst (randomly)
                    c = random.choice(self.catalyst_set)
                    cx= bin(c)
                    cx= cx[2:]
                    react['c'] = cx.rjust(self.strng_size, '.')
                    self.ligand_reactions.append(react)
                    # include also the opposite
                    # reaction
                    if x1 != x2:
                        y2 = x2 + x1
                        # convert to decimal
                        pr2= self.convert_binstr_to_dec(y2)
                        if pr2 != pr:
                            # set new dictionary
                            react = {}
                            react['r1'] = x2.rjust(self.strng_size, '.')
                            react['r2'] = x1.rjust(self.strng_size, '.')
                            react['p']  = y2.rjust(self.strng_size, '.')
                            # select catalyst (randomly)
                            c2 = random.choice(self.catalyst_set)
                            cx= bin(c2)
                            cx= cx[2:]
                            react['c'] = cx.rjust(self.strng_size, '.')
                            self.ligand_reactions.append(react)
        if log.level <= logging.DEBUG:
            log.info("\n")
            for r in self.ligand_reactions:
                log.debug("\t r1 : " + r['r1'] + " |\t r2 : " + r['r2'] + " |\t c : " + r['c'] + " |\t p : " + r['p'])
            log.info("\t " + p.sep)
        # for cleavage
        # 1) 'r'  : reactant
        # 2) 'c'  : catalyst
        # 3) 'p1' : product 1
        # 4) 'p2' : product 2
        tmp_list = []
        for r in range(1, self.size_X+1):
            x = bin(r)
            x = x[2:]
            n = len(x)
            if n > 1:
                for i in range(1, n):
                    x1 = x[:i]
                    p1 = self.convert_binstr_to_dec(x1)
                    x2 = x[i:]
                    p2 = self.convert_binstr_to_dec(x2)
                    if p1 > 0 and p2 > 0:
                        k=0
                        while x1[k] == '0':
                            k+=1
                        x1 = x1[k:]
                        k=0
                        while x2[k] == '0':
                            k+=1
                        x2 = x2[k:]
                        tmp_list.append([x, x1, x2])
        # screen list of eqv. reactions
        to_remove = []
        i = 0
        while i < len(tmp_list):
            [x, x1, x2] = tmp_list[i]
            for j in range(i+1, len(tmp_list)):
                [y, y1, y2] = tmp_list[j]
                if x == y and x1 == y2 and x2 == y1:
                    to_remove.append(j)
                    i = j+1
                    break
            if j == len(tmp_list)-1:
                i += 1
        react_lst = []
        for i in range(len(tmp_list)):
            if i in to_remove:
                pass
            else:
                react_lst.append(tmp_list[i])
        # build dict.
        for r in range(len(react_lst)):
            [x, x1, x2] = react_lst[r]
            # set new dictionary
            react = {}
            react['r'] = x.rjust(self.strng_size, '.')
            react['p1']= x1.rjust(self.strng_size, '.')
            react['p2']= x2.rjust(self.strng_size, '.')
            # select catalyst (randomly)
            c = random.choice(self.catalyst_set)
            cx= bin(c)
            cx= cx[2:]
            react['c'] = cx.rjust(self.strng_size, '.')
            self.cleavage_reactions.append(react)
        if log.level <= logging.DEBUG:
            log.info("\n")
            for r in self.cleavage_reactions:
                log.debug("\t r : " + r['r'] + " |\t c : " + r['c'] + " |\t p1 : " + r['p1'] + " |\t p2 : " + r['p2'])
            log.info("\t " + p.sep)