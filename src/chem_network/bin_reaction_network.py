import math
import random
import numpy as np
import logging
import os
from src.utilities.logging_module import log
from src.input_data.read_input import p
from src.utilities.graph_class import graph_obj
from src.molecules_dyn.gillespie_algo import chemical_kinetics_solver
from src.catalysts.catalysts_set import build_catalysts_list
from src.utilities.plot_catal_distr import plot_ACFS_hist_distr
from src.chem_network.network_fitness import network_fitness
from src.chem_network.abstract_reaction_network import ReactionNetwork
from src.reactions.reaction_class import LigationReaction, CleavageReaction

#  class describing the
#  reaction network -> we use a binary polymer model


class BinaryReactionNetwork(ReactionNetwork):
    """
    Binary polymer reaction network without catalyst assignment.

    This class only builds the reaction topology.
    Catalyst assignment is handled later by the protocell.
    """
    def __init__(self):
        super().__init__(
            network_type="binary_polymer"
        )
        self.reaction_list = []
    # --------------------------------------------------
    # required interface
    # --------------------------------------------------
    def build_network(self, X_set):
        self._clear_reactions()
        self._build_species_maps(X_set)
        self._build_ligation_reactions()
        self._build_cleavage_reactions()
        self._validate_network()
    # --------------------------------------------------
    # validation
    # --------------------------------------------------
    def _validate_network(self):
        """
        Validate internal consistency of the reaction network object.
        Raise ValueError with a clear message on failure.
        """
        # basic sizes
        if self.size_X <= 0:
            log.error("size_X must be > 0")
        if not hasattr(self, "species_set") or len(self.species_set) == 0:
            log.error("species_set is empty")
        # species maps
        if not hasattr(self, "X_set") or len(self.X_set) != self.size_X:
            log.error("X_set size is inconsistent with size_X")
        if not hasattr(self, "x_to_id") or not hasattr(self, "id_to_x"):
            log.error("species maps x_to_id / id_to_x are missing")
        if len(self.x_to_id) != self.size_X:
            log.error("x_to_id size is inconsistent with size_X")
        if len(self.id_to_x) != self.size_X:
            log.error("id_to_x size is inconsistent with size_X")
        # species ids must match convention 1..size_X
        expected_species = list(range(1, self.size_X + 1))
        if list(self.species_set) != expected_species:
            log.error(
                f"species_set is inconsistent: expected {expected_species}, got {self.species_set}"
            )
        # reverse consistency of maps
        for sid in self.species_set:
            if sid not in self.id_to_x:
                log.error(f"id {sid} missing in id_to_x")
            mol = self.id_to_x[sid]
            if mol not in self.x_to_id:
                log.error(f"molecule for id {sid} missing in x_to_id")
            if self.x_to_id[mol] != sid:
                log.error(
                    f"inconsistent species maps for id {sid}: x_to_id[id_to_x[{sid}]] != {sid}"
                )
        # max string size consistency
        max_len = max(len(mol) for mol in self.X_set)
        if self.max_strng_size != max_len:
            log.error(
                f"max_strng_size inconsistent: expected {max_len}, got {self.max_strng_size}"
            )
        # sequence lookup map if present
        if hasattr(self, "seq_to_id"):
            if len(self.seq_to_id) != self.size_X:
                log.error("seq_to_id size is inconsistent with size_X")
            for sid, mol in self.id_to_x.items():
                seq = mol.show_sequence()
                if seq not in self.seq_to_id:
                    log.error(f"sequence '{seq}' missing in seq_to_id")
                if self.seq_to_id[seq] != sid:
                    log.error(
                        f"inconsistent seq_to_id for sequence '{seq}': expected id {sid}, got {self.seq_to_id[seq]}"
                    )
        # reaction list if present
        if hasattr(self, "reaction_list"):
            seen = set()
            for rxn in self.reaction_list:
                if not hasattr(rxn, "reaction_id"):
                    log.error("reaction without reaction_id found")
                if not hasattr(rxn, "reaction_type"):
                    log.error(f"reaction {rxn.reaction_id} missing reaction_type")
                if rxn.reaction_type not in {"ligation", "cleavage"}:
                    log.error(
                        f"reaction {rxn.reaction_id} has invalid reaction_type '{rxn.reaction_type}'"
                    )
                if not hasattr(rxn, "reactants") or not hasattr(rxn, "products"):
                    log.error(
                        f"reaction {rxn.reaction_id} missing reactants/products"
                    )
                # check species ids in reactions
                for sid in list(rxn.reactants) + list(rxn.products):
                    if sid not in self.id_to_x:
                        log.error(
                            f"reaction {rxn.reaction_id} references unknown species id {sid}"
                        )
                # uniqueness check
                key = (rxn.reaction_type, tuple(rxn.reactants), tuple(rxn.products))
                if key in seen:
                    log.error(
                        f"duplicate reaction found: type={rxn.reaction_type}, "
                        f"reactants={rxn.reactants}, products={rxn.products}"
                    )
                seen.add(key)
        log.info("\n")
        log.info("\t NETWORK VALIDATION TEST PASSED")
        log.info("\t " + p.sep)
    # --------------------------------------------------
    # species set
    # --------------------------------------------------
    def _build_species_maps(self, X_set):
        """
        Build mappings:
            molecule string -> molecule id
            molecule id     -> molecule string

        Molecule ids are 1-based to match your current convention.
        """
        self.X_set = list(X_set)
        self.x_to_id = {}
        self.id_to_x = {}
        self.seq_to_id = {}
        # run over X_set
        for i, mol in enumerate(self.X_set, start=1):
            self.x_to_id[mol] = i
            self.id_to_x[i] = mol
            self.seq_to_id[mol.show_sequence()] = i
        # species set
        self.species_set = list(range(1, len(self.X_set) + 1))
        self.size_X = len(X_set)
        self.max_strng_size = max(len(mol) for mol in X_set)
    # --------------------------------------------------
    # reaction builders
    # --------------------------------------------------
    def _build_ligation_reactions(self):
        reaction_id = 0
        stored = set()
        for r1 in self.species_set:
            x1 = self.id_to_x[r1]
            n1 = len(x1)
            for r2 in self.species_set:
                x2 = self.id_to_x[r2]
                n2 = len(x2)
                if n1 + n2 > self.max_strng_size:
                    continue
                y = x1.ligate(x2)
                y_seq = y.show_sequence()
                # 1) convert back to species id
                if y_seq not in self.seq_to_id:
                    continue
                p_id = self.seq_to_id[y_seq]
                # ligation key
                rl_key = (r1, r2, p_id)
                if rl_key in stored:
                    continue
                stored.add(rl_key)
                react = LigationReaction(
                    r1=r1, 
                    r2=r2, 
                    p=p_id, 
                    reaction_id=reaction_id
                )
                self.add_reaction(react)
                reaction_id += 1
    # cleavage reactions
    def _build_cleavage_reactions(self):
        if len(self.reaction_list) == 0:
            reaction_id = 0
        else:
            reaction_id = self.reaction_list[-1].reaction_id + 1
        stored = set()
        for r in self.species_set:
            x = self.id_to_x[r]
            n = len(x)
            if n <= 1:
                continue
            x_seq = x.show_sequence()
            for i in range(1, n):
                # split the sequence
                s1 = x_seq[:i]
                s2 = x_seq[i:]
                if len(s1) == 0 or len(s2) == 0:
                    continue
                # map products to id
                if s1 not in self.seq_to_id or s2 not in self.seq_to_id:
                    continue
                p1_id = self.seq_to_id[s1]
                p2_id = self.seq_to_id[s2]
                # cleavage key
                rc_key = (r, p1_id, p2_id)
                if rc_key in stored:
                    continue
                stored.add(rc_key)
                react = CleavageReaction(
                    r=r,
                    p1=p1_id, 
                    p2=p2_id, 
                    reaction_id=reaction_id
                )
                self.add_reaction(react)
                reaction_id += 1
    # --------------------------------------------------
    # utilities
    # --------------------------------------------------
    def convert_binstr_to_dec(self, y):
        return int(y.zfill(self.strng_size), 2)
    def convert_dec_to_binstr(self, x):
        return bin(x)[2:]
    def fmt_str(self, x):
        return x.rjust(self.strng_size, ".")
    # --------------------------------------------------
    # convenience
    # --------------------------------------------------
    def get_reaction_by_id(self, reaction_id):
        for rxn in self.reaction_list:
            if rxn.reaction_id == reaction_id:
                return rxn
        log.error(f"Reaction id {reaction_id} not found")
    def summary(self):
        n_lig = sum(1 for r in self.reaction_list if r.reaction_type == "ligation")
        n_clv = sum(1 for r in self.reaction_list if r.reaction_type == "cleavage")
        log.info(f"\t network_type: {self.network_type}")
        log.info(f"\t max strng_size: {self.max_strng_size}")
        log.info(f"\t size_X: {self.size_X}")
        log.info(f"\t n_species: {len(self.species_set)}")
        log.info(f"\t n_reactions: {len(self.reaction_list)}")
        log.info(f"\t n_ligations: {n_lig}")
        log.info(f"\t n_cleavages: {n_clv}")
        log.info("\t " + p.sep)















class _BinaryReactionNetwork:
    """
    Binary polymer model reaction network class
    """
    def __init__(self):
        # max. num. reactions
        self.num_reactions = None
        # reactions set
        self.ligand_reactions = None
        self.cleavage_reactions = None
    def set_binary_polymer_model(self, catalyst_distr):
        # n. food set bits
        self.n_F_bits = math.log2(self.size_F)
        log.info("\t max. string size : " + str(self.strng_size))
        # build the catalysts set (C)
        self.define_catalysts_set(catalyst_distr)
        # build food set
        self.build_food_set()
        # build set of reactions
        self.build_reactions_set()
    '''
    def set_binary_polymer_from_genome(self, catalyst_distr):
        # n. food set bits
        self.n_F_bits = math.log2(self.size_F)
        log.info("\t max. string size : " + str(self.strng_size))
        # build catalysts set C
        self.define_catalysts_set(catalyst_distr)
        # set catalysts
        catal_set_react = self.genome_to_catalysts()
        # builid food set
        self.build_food_set()
        # build reactions
        self.build_reactions_set_from_catalysts(catal_set_react)
    '''
    # catalysts set
    def define_catalysts_set(self, catalyst_distr):
        self.catalyst_distr = catalyst_distr
        log.info("\n")
        log.info("\t C = " + str(self.catalyst_distr))
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
                    # select catalyst list from distribution
                    if self.catalyst_distr is None:
                        c = random.choice(self.catalyst_set)
                        react['c_int'] = c
                        cx= bin(c)
                        cx= cx[2:]
                        react['c'] = cx.rjust(self.strng_size, '.')
                    else:
                        react['cList'] = build_catalysts_list(self.catalyst_distr, p.nconfig)
                        if log.level <= logging.INFO:
                            output_file = p.working_dir + "/ACF_HIST_catal_distr" + str(r1) + "-" + str(r2) + "-" + str(self.typ_index) + "-" + str(self.net_index) + ".pdf"
                            plot_ACFS_hist_distr(react['cList'], self.size_X, output_file)
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
                            if self.catalyst_distr is None:
                                c2 = random.choice(self.catalyst_set)
                                react['c_int'] = c2
                                cx= bin(c2)
                                cx= cx[2:]
                                react['c'] = cx.rjust(self.strng_size, '.')
                            else:
                                react['cList'] = build_catalysts_list(self.catalyst_distr, p.nconfig)
                            self.ligand_reactions.append(react)
        if log.level <= logging.DEBUG:
            log.info("\n")
            for r in self.ligand_reactions:
                if self.catalyst_distr is None:
                    log.debug("\t r1 : " + r['r1'] + " |\t r2 : " + r['r2'] + " |\t c : " + r['c'] + " |\t p : " + r['p'])
                else:
                    log.debug("\t r1 : " + r['r1'] + " |\t r2 : " + r['r2'] + " |\t c : " + r['cList'] + " |\t p : " + r['p'])
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
            # select catalyst distribution
            if self.catalyst_distr is None:
                c = random.choice(self.catalyst_set)
                react['c_int'] = c
                cx= bin(c)
                cx= cx[2:]
                react['c'] = cx.rjust(self.strng_size, '.')
            else:
                react['cList'] = build_catalysts_list(self.catalyst_distr, p.nconfig)
            self.cleavage_reactions.append(react)
        if log.level <= logging.DEBUG:
            log.info("\n")
            for r in self.cleavage_reactions:
                if self.catalyst_distr is None:
                    log.debug("\t r : " + r['r'] + " |\t c : " + r['c'] + " |\t p1 : " + r['p1'] + " |\t p2 : " + r['p2'])
                else:
                    log.debug("\t r : " + r['r'] + " |\t c : " + r['cList'] + " |\t p1 : " + r['p1'] + " |\t p2 : " + r['p2'])
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
    #
    # define the reaction kinetic
    # model
    def set_fitness_from_chemical_kinetics(self, nkin_simul, molecules_fitness, fitness_p, max_fitness):
        # first set the solver
        kinetic_solver = chemical_kinetics_solver(nkin_simul)
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
        # run over catalysts list
        List_state_t = [None]*p.nconfig
        List_times = [None]*p.nconfig
        for ic in range(p.nconfig):
            log.info("\t configuration: " + str(ic) + " / " + str(p.nconfig))
            kinetic_solver.set_propensity(reaction_set, self.size_X, ic)
            t_grid, avg_states_t = kinetic_solver.solve()
            List_times[ic] = t_grid
            List_state_t[ic] = avg_states_t
        t_grid, avg_states_t = kinetic_solver.average_trajectories(List_times, List_state_t)
        # define network fitness
        self.netw_fitness = network_fitness(kinetic_solver.X_set, 
                                            kinetic_solver.X_mass, 
                                            molecules_fitness, 
                                            max_fitness)
        if log.level == logging.INFO:
            self.netw_fitness.show_fitness_distr()
        fitness = self.netw_fitness.set_fitness(t_grid, avg_states_t, fitness_p)
        log.info("\t fitness value: " + str(fitness))
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