import math
import random
import numpy as np
import logging
from lifeorig.logging_module import log
from lifeorig.read_input import p
from lifeorig.graph_class import graph_obj
from lifeorig.gillespie_algo import chemical_kinetics_solver
import gravis as gv
import os
#  class describing the
#  reaction network -> we use a binary polymer model
class reaction_net_class:
    def __init__(self, size_bpol, size_F, size_C):
        # max. size string
        # polymer model
        self.strng_size = size_bpol
        # food set size
        self.size_F = size_F
        self.food_set = []
        # catalysts size
        self.size_C = size_C
        self.catalyst_set = []
        # size of molecules set
        self.size_X = 2 ** size_bpol - 1
        # reactions set
        self.ligand_reactions = []
        self.cleavage_reactions = []
        # fitness
        self.fitness = 0.
    def set_binary_polymer_model(self, catalyst_set):
        # n. food set bits
        self.n_F_bits = math.log2(self.size_F)
        log.info("\t max. string size : " + str(self.strng_size))
        # build the catalysts set (C)
        self.define_catalysts_set(catalyst_set)
        # builid food set
        self.build_food_set()
        # build set of reactions
        self.build_reactions_set()
    def set_binary_polymer_from_genome(self, catalyst_set):
        # n. food set bits
        self.n_F_bits = math.log2(self.size_F)
        log.info("\t max. string size : " + str(self.strng_size))
        # build catalysts set C
        self.define_catalysts_set(catalyst_set)
        # set catalysts
        catal_set_react = self.genome_to_catalysts()
        # builid food set
        self.build_food_set()
        # build reactions
        self.build_reactions_set_from_catalysts(catal_set_react)
    # catalysts set
    def define_catalysts_set(self, catalyst_set):
        self.catalyst_set = catalyst_set
        log.info("\n")
        log.info("\t C = " + str(self.catalyst_set))
        log.info("\t " + p.sep)
    # build food set
    def build_food_set(self):
        for i in range(1, self.size_F+1):
            self.food_set.append(i)
        log.info("\n")
        log.info("\t F = " + str(self.food_set))
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
                    react['r1_int'] = self.convert_binstr_to_dec(x1)
                    react['r2'] = x2.rjust(self.strng_size, '.')
                    react['r2_int'] = self.convert_binstr_to_dec(x2)
                    react['p']  = y.rjust(self.strng_size, '.')
                    react['p_int'] = self.convert_binstr_to_dec(y)
                    # select catalyst (randomly)
                    c = random.choice(self.catalyst_set)
                    react['c_int'] = c
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
                            react['r1_int'] = self.convert_binstr_to_dec(x2)
                            react['r2'] = x1.rjust(self.strng_size, '.')
                            react['r2_int'] = self.convert_binstr_to_dec(x1)
                            react['p']  = y2.rjust(self.strng_size, '.')
                            react['p_int'] = self.convert_binstr_to_dec(y2)
                            # select catalyst (randomly)
                            c2 = random.choice(self.catalyst_set)
                            react['c_int'] = c2
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
            react['r_int'] = self.convert_binstr_to_dec(x)
            react['p1']= x1.rjust(self.strng_size, '.')
            react['p1_int']= self.convert_binstr_to_dec(x1)
            react['p2']= x2.rjust(self.strng_size, '.')
            react['p2_int']= self.convert_binstr_to_dec(x2)
            # select catalyst (randomly)
            c = random.choice(self.catalyst_set)
            react['c_int'] = c
            cx= bin(c)
            cx= cx[2:]
            react['c'] = cx.rjust(self.strng_size, '.')
            self.cleavage_reactions.append(react)
        if log.level <= logging.DEBUG:
            log.info("\n")
            for r in self.cleavage_reactions:
                log.debug("\t r : " + r['r'] + " |\t c : " + r['c'] + " |\t p1 : " + r['p1'] + " |\t p2 : " + r['p2'])
            log.info("\t " + p.sep)
    #
    # set network genome
    def set_network_genome(self):
        self.genome = ""
        # first ligand reactions
        for r in self.ligand_reactions:
            c = r['c_int']
            cx= bin(c)
            cx= cx[2:]
            self.genome += cx.zfill(self.strng_size)
        for r in self.cleavage_reactions:
            c = r['c_int']
            cx= bin(c)
            cx= cx[2:]
            self.genome += cx.zfill(self.strng_size)
        #
        if log.level <= logging.INFO:
            log.info("\n")
            log.info("\t network genome : " + self.genome)
            log.info("\n")
            log.info("\t " + p.sep)
    #
    # build reaction set given
    # catalyst set (input)
    def build_reactions_set_from_catalysts(self, catalysts):
        # list of dictionaries with keys
        # for ligation
        # 1) 'r1' : first reactant
        # 2) 'r2' : second reactant
        # 3) 'c'  : catalyst
        # 4) 'p'  : product
        c_i = 0
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
                    react['r1_int'] = self.convert_binstr_to_dec(x1)
                    react['r2'] = x2.rjust(self.strng_size, '.')
                    react['r2_int'] = self.convert_binstr_to_dec(x2)
                    react['p']  = y.rjust(self.strng_size, '.')
                    react['p_int'] = self.convert_binstr_to_dec(y)
                    # select catalyst
                    c = catalysts[c_i]
                    c_i += 1
                    react['c_int'] = c
                    cx = bin(c)
                    cx = cx[2:]
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
                            react['r1_int'] = self.convert_binstr_to_dec(x2)
                            react['r2'] = x1.rjust(self.strng_size, '.')
                            react['r2_int'] = self.convert_binstr_to_dec(x1)
                            react['p']  = y2.rjust(self.strng_size, '.')
                            react['p_int'] = self.convert_binstr_to_dec(y2)
                            # select catalyst
                            c2 = catalysts[c_i]
                            c_i += 1
                            react['c_int'] = c2
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
            react['r_int'] = self.convert_binstr_to_dec(x)
            react['p1']= x1.rjust(self.strng_size, '.')
            react['p1_int']= self.convert_binstr_to_dec(x1)
            react['p2']= x2.rjust(self.strng_size, '.')
            react['p2_int']= self.convert_binstr_to_dec(x2)
            # select catalyst
            c = catalysts[c_i]
            c_i += 1
            react['c_int'] = c
            cx= bin(c)
            cx= cx[2:]
            react['c'] = cx.rjust(self.strng_size, '.')
            self.cleavage_reactions.append(react)
        if log.level <= logging.DEBUG:
            log.info("\n")
            for r in self.cleavage_reactions:
                log.debug("\t r : " + r['r'] + " |\t c : " + r['c'] + " |\t p1 : " + r['p1'] + " |\t p2 : " + r['p2'])
            log.info("\t " + p.sep)
    def set_genome_input(self, genome):
        self.genome = ""
        c_i = 0
        while(c_i < len(genome)):
            cat_bin = genome[c_i:c_i+self.strng_size]
            c_i += self.strng_size
            # compute catalyst
            r = 0
            j = 0
            for c in cat_bin:
                r += int(c)* 2**(self.strng_size-1-j)
                j += 1
            if r == 0:
                log.warning("\t catalyst -> 0")
                log.warning("\t set to 1")
                self.genome += "0001"
            else:
                self.genome += cat_bin
    #
    # convert network genome to catalysts
    # list
    def genome_to_catalysts(self):
        c_i = 0
        catalysts_lst = []
        while (c_i < len(self.genome)):
            cat_bin = self.genome[c_i:c_i+self.strng_size]
            c_i += self.strng_size
            r = 0
            j = 0
            for c in cat_bin:
                r += int(c)* 2**(self.strng_size-1-j)
                j += 1
            catalysts_lst.append(r)
        return catalysts_lst
    def compute_fitness(self, kinetic_solver):
        avg_state_t = kinetic_solver.avg_state_t
        # fitness = \sum_target molecules ni / N
        T = -1
        target_molecules = p.target_molecules
        # molecular masses
        X_mass = kinetic_solver.X_mass
        # n. target molecules
        nt = len(target_molecules)
        target_mass = 0
        for i in range(nt):
            ml = target_molecules[i]
            target_mass += avg_state_t[T,ml] * X_mass[ml]
        # total mass
        total_mass = 0
        for i in range(self.size_X):
            total_mass += avg_state_t[T,i] * X_mass[i]
        self.fitness = target_mass / total_mass
        log.info("\n")
        log.info("\t " + p.sep)
        log.info("\t network fitness : " + str(self.fitness))
        log.info("\t " + p.sep)
    #
    # define the reaction kinetic
    # model
    def set_chemical_kinetics_solver(self):
        # first set the solver
        kinetic_solver = chemical_kinetics_solver()
        # reaction set full list
        reaction_set = []
        for r in self.ligand_reactions:
            reaction_set.append(r)
        for r in self.cleavage_reactions:
            reaction_set.append(r)
        # set stoichiometry
        # list
        kinetic_solver.build_X_set(self.size_X)
        kinetic_solver.set_initial_population(self.food_set)
        kinetic_solver.set_stoichiometry(reaction_set, self.size_X)
        kinetic_solver.set_propensity(reaction_set, self.size_X)
        kinetic_solver.solve()
        # molecules to display
        target_molecules = p.target_molecules
        if log.level == logging.DEBUG:
            kinetic_solver.show(target_molecules)
        # set fitness of chemical
        # network
        self.compute_fitness(kinetic_solver)
    #
    # find ACF subset
    # this subroutine find RAF subset if present in the network
    # RAF subset has 2 properties
    # 1- must be autocatalytic
    # 2- food generated
    def find_ACF_subset(self):
        Nl = len(self.ligand_reactions)
        Nc = len(self.cleavage_reactions)
        R2 = list(np.arange(Nl+Nc))
        R  = []
        F = np.arange(1, self.size_F+1)
        # reaction set - full list
        reaction_set = []
        for r in self.ligand_reactions:
            reaction_set.append(r)
        for r in self.cleavage_reactions:
            reaction_set.append(r)
        # find final reaction set
        while R2 != R:
            R = R2
            R2 = []
            # compute Cl_R(F)
            Cl_R = self.compute_closure_set(R, F, reaction_set)
            # find new reaction set
            for ri in R:
                r = reaction_set[ri]
                if 'r1_int' in r and 'r2_int' in r:
                    r1 = r['r1_int']
                    r2 = r['r2_int']
                    c  = r['c_int']
                    if r1 in Cl_R and r2 in Cl_R and c in Cl_R:
                        R2.append(ri)
                elif 'r_int' in r:
                    r1 = r['r_int']
                    c  = r['c_int']
                    if r1 in Cl_R and c in Cl_R:
                        R2.append(ri)
            #print(R, R2)
        # set up ACF set
        self.ACF_set = []
        for ri in R:
            r = reaction_set[ri]
            self.ACF_set.append(r)
    # routine to compute
    # the closure set of F
    def compute_closure_set(self, R, F, reaction_set):
        # input 1 : R - set of reactions
        # input 2 : X - molecules set
        X = set(F)
        Y = set()
        while Y != X:
            Y = X.copy()
            for r in reaction_set:
                if reaction_set.index(r) in R:
                    # ligand reactions
                    if 'r1_int' in r and 'r2_int' in r:
                        r1 = r['r1_int']
                        r2 = r['r2_int']
                        c  = r['c_int']
                        if r1 in X and r2 in X and c in X:
                            X.add(r['p_int'])
                    # cleavage reactions
                    elif 'r_int' in r:
                        r1 = r['r_int']
                        c  = r['c_int']
                        if r1 in X and c in X:
                            X.add(r['p1_int'])
                            X.add(r['p2_int'])
        return Y
    #
    # show the reaction network
    def show_network_test(self):
        graph1 = {
            'graph':{
                'directed': True,
                'metadata': {
                    'arrow_size': 5,
                    'background_color': 'black',
                    'edge_size': 3,
                    'edge_label_size': 14,
                    'edge_label_color': 'white',
                    'node_size': 15,
                    'node_color': 'white',
                },
                'nodes': {
                    1: {'metadata': {'shape': 'rectangle', 'y': 200}},
                    2: {},
                    3: {},
                    4: {'metadata': {'shape': 'rectangle', 'y': 200}},
                    5: {'metadata': {'shape': 'hexagon', 'y': 0}},
                },
                'edges': [
                    {'source': 1, 'target': 2, 'metadata': {'color': '#d73027', 'de': 'Das',   'en': 'This'}},
                    {'source': 2, 'target': 3, 'metadata': {'color': '#f46d43', 'de': 'ist',   'en': 'is'}},
                    {'source': 3, 'target': 1, 'metadata': {'color': '#fdae61', 'de': 'das',   'en': 'the'}},
                    {'source': 1, 'target': 4, 'metadata': {'color': '#fee08b', 'de': 'Haus',  'en': 'house'}},
                    {'source': 4, 'target': 3, 'metadata': {'color': '#d9ef8b', 'de': 'vom',   'en': 'of'}},
                    {'source': 3, 'target': 5, 'metadata': {'color': '#a6d96a', 'de': 'Ni-.',  'en': 'San-'}},
                    {'source': 5, 'target': 2, 'metadata': {'color': '#66bd63', 'de': 'ko-',   'en': 'ta'}},
                    {'source': 2, 'target': 4, 'metadata': {'color': '#1a9850', 'de': 'laus.', 'en': 'Claus.'}},
                ],
            }
        }
        fig = gv.three(graph1, show_node_label=False, show_edge_label=True, edge_label_data_source='en')
        isExist = os.path.isfile(p.working_dir+'/g1.html')
        if not isExist:
            fig.export_html(p.working_dir+'/g1.html')
        #
        graph_gjgf = {
            'graph': {
                'directed': True,
                'metadata': {
                    'node_label_size': 14,
                    'node_label_color': 'green',
                    'edge_label_size': 10,
                    'edge_label_color': 'blue',
                },
                'nodes': [
                    {'id': '0', 'label': 'first node', 'metadata': {
                    'color': 'red',
                    'size': 15,
                    'shape': 'rectangle',
                    'opacity': 0.7,
                    'label_color': 'red',
                    'label_size': 20,
                    'border_color': 'black',
                    'border_size': 3,
                    }},
                    {'id': '1'},
                    {'id': '2'},
                    {'id': '3', 'metadata': {
                        'color': 'green',
                        'size': 15,
                        'shape': 'hexagon',
                        'opacity': 0.7,
                        'label_color': 'green',
                        'label_size': 10,
                        'border_color': 'blue',
                        'border_size': 3,
                    }},
                    {'id': '4'},
                    {'id': '5'},
                    {'id': '6', 'label': 'last node'},
                ],
                'edges': [
                    {'source': 0, 'target': 1},
                    {'source': 1, 'target': 2, 'label': 'e2'},
                    {'source': 2, 'target': 3},
                    {'source': 3, 'target': 4},
                    {'source': 4, 'target': 5, 'label': 'e5', 'metadata': {
                        'color': 'orange',
                        'label_color': 'gray',
                        'label_size': 14,
                        'size': 4.0,
                    }},
                    {'source': 5, 'target': 6},
                    {'source': 6, 'target': 2, 'label': 'e2'},
                ]
            }
        }
        fig = gv.three(graph_gjgf, graph_height=200,
            node_label_data_source='label',
            show_edge_label=True, edge_label_data_source='label')
        isExist = os.path.isfile(p.working_dir+'/g2.html')
        if not isExist:
            fig.export_html(p.working_dir+'/g2.html')
    #
    # show full network
    def show_network(self, file_name):
        graph = graph_obj()
        graph.set_metadata(14, 'green', 10, 'blue')
        # first add the food set nodes
        for i in self.food_set:
            label = 'f'+str(i)
            graph.add_node(i-1, label, 'blue', 10, 'circle', 'black', 20, 'black', 3)
        # add remaining chemicals
        for i in range(1, self.size_X+1):
            if i not in self.food_set:
                label = 'p'+str(i)
                graph.add_node(i-1, label, 'black', 10, 'circle', 'black', 20, 'black', 3)
        # add edges
        react_index = self.size_X+1
        for r in self.ligand_reactions:
            label = 'r'+str(react_index-self.size_X)
            if r not in self.ACF_set:
                graph.add_node(react_index-1, label, 'white', 10, 'rectangle', 'black', 20, 'black', 3)
            else:
                graph.add_node(react_index-1, label, 'white', 10, 'rectangle', 'green', 20, 'green', 3)
            graph.add_edge(r['r1_int']-1, react_index-1, 'black')
            if r['r2_int'] != r['r1_int']:
                graph.add_edge(r['r2_int']-1, react_index-1, 'black')
            graph.add_edge(react_index-1, r['p_int']-1, 'black')
            graph.add_edge(int(r['c_int'])-1, react_index-1, 'red')
            react_index += 1
        # add cleavage reactions
        for r in self.cleavage_reactions:
            label = 'r'+str(react_index-self.size_X)
            if r not in self.ACF_set:
                graph.add_node(react_index-1, label, 'white', 10, 'rectangle', 'black', 20, 'black', 3)
            else:
                graph.add_node(react_index-1, label, 'white', 10, 'rectangle', 'green', 20, 'green', 3)
            graph.add_node(react_index-1, label, 'white', 10, 'rectangle', 'black', 20, 'black', 3)
            graph.add_edge(r['r_int']-1, react_index-1, 'black')
            graph.add_edge(react_index-1, r['p1_int']-1, 'black')
            if r['p2_int'] != r['p1_int']:
                graph.add_edge(react_index-1, r['p2_int']-1, 'black')
            graph.add_edge(int(r['c_int'])-1, react_index-1, 'red')
            react_index += 1
        # prepare figure
        fig = gv.d3(graph.graph_gjgf, graph_height=300,
            node_label_data_source='label',
            show_edge_label=True, edge_label_data_source='label',
            edge_curvature=0.4, zoom_factor=2.5)
        isExist = os.path.isfile(file_name)
        if not isExist:
            fig.export_html(file_name)