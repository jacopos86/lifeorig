import numpy as np
import math
from src.read_input import p
from src.logging_module import log
from src.reaction_network import reaction_net_class
# This class defines the fitness
# function
class fitness_distr():
    def __init__(self, n):
        # set fitness array
        self.fitness = np.zeros(n)
        self.fitness_oft = None
        # shape
        self.size = n
    def set_fitness_distr(self, ACF_set):
        # set up distribution
        if p.fitness_eval == "random":
            self.set_random_initial_fitness(p.seed, p.max_fitness)
        elif p.fitness_eval == "compute":
            i = 0
            for ACF in ACF_set:
                self.fitness[i] = ACF.fitness
                i += 1
        elif p.fitness_eval == "read":
            # open data file
            file_name = p.working_dir + "/ACF_data.txt"
            self.extract_fitness_from_file(file_name, ACF_set)
        else:
            log.error("fitness_eval not recognized")
    def set_random_initial_fitness(self, s, M):
        np.random.seed(s)
        # set [0, M) random distribution
        self.fitness = np.random.rand(self.size)
        self.fitness[:] = self.fitness[:] * M
    def extract_fitness_from_file(self, file_name, ACF_set):
        # open file
        f = open(file_name, 'r')
        lines = f.readlines()
        assert len(lines) == len(ACF_set)
        # extract data
        i = 0
        for line in lines:
            line = line.strip().split()
            # set genome
            ACF_set[i].genome = ''
            ACF_set[i].genome = line[0]
            # fitness
            ACF_set[i].fitness= float(line[2])
            self.fitness[i] = ACF_set[i].fitness
            i += 1
    def set_constant_fitness_over_time(self, nt):
        self.fitness_oft = np.zeros((self.size,nt))
        for t in range(nt):
            self.fitness_oft[:,t] = self.fitness[:]
    def show_fitness_distr(self):
        log.info("\t " + p.sep)
        log.info("\n")
        log.info("\t fitness distr : ")
        log.info("\n")
        line = ""
        for j in range(self.size):
            line += " {0:.3f}".format(self.fitness[j])
        log.info("\t " + line)
        log.info("\n")
        log.info("\t " + p.sep)
        
#
#  fitness distribution
#  for evolutionary game dynamics
#   f_i = \sum_j a_ij x_j

class fitness_distr_game_dyn():
    def __init__(self, n):
        # payoff matrix
        self.a_ij = None
        self.aij_oft = None
        # shape
        self.size = n
    # set fitness distr.
    def set_fitness_distr(self, ACF_set, out_file):
        # set distribution
        if p.fitness_eval == "compute":
            self.set_payoff_matrix(ACF_set, out_file)
        elif p.fitness_eval == "read":
            # open file
            inp_file = out_file
            self.extract_fitness_from_file(inp_file, ACF_set)
        else:
            log.error("fitness_eval not recognized")
    # define payoff matrix
    def set_payoff_matrix(self, ACF_distr, out_file):
        # initialize payoff
        self.a_ij = np.zeros((self.size, self.size))
        # define new reaction network
        # new network = ACFd[i] + ACFd[j]
        # join two sets of reactions
        for i in range(self.size):
            ACFS_i = ACF_distr[i]
            self.a_ij[i,i] = ACFS_i.fitness
            for j in range(i+1, self.size):
                ACFS_j = ACF_distr[j]
                # build total
                # reaction set
                #print(ACFS_i.genome_to_catalysts())
                ACFS_ij = reaction_net_class(p.bpol_strng_size, p.size_F, p.size_C)
                # catalyst list
                catalyst_set = ACFS_i.catalyst_set + ACFS_j.catalyst_set
                catalyst_set = list(set(catalyst_set))
                # n. food set bits
                ACFS_ij.n_F_bits = math.log2(ACFS_ij.size_F)
                log.info("\t max. string size : " + str(ACFS_ij.strng_size))
                # build the catalysts set (C)
                ACFS_ij.define_catalysts_set(catalyst_set)
                # builid food set
                ACFS_ij.build_food_set()
                # build reaction set
                for r_i in ACFS_i.ligand_reactions:
                    ACFS_ij.ligand_reactions.append(r_i)
                for r_j in ACFS_j.ligand_reactions:
                    for r_ij in ACFS_ij.ligand_reactions:
                        if r_j['r1_int'] == r_ij['r1_int'] and r_j['r2_int'] == r_ij['r2_int'] and r_j['p_int'] == r_ij['p_int']:
                            if r_j['c_int'] != r_ij['c_int']:
                                ACFS_ij.ligand_reactions.append(r_j)
                            else:
                                break
                for r_i in ACFS_i.cleavage_reactions:
                    ACFS_ij.cleavage_reactions.append(r_i)
                for r_j in ACFS_j.cleavage_reactions:
                    for r_ij in ACFS_ij.cleavage_reactions:
                        if r_j['r_int'] == r_ij['r_int'] and r_j['p1_int'] == r_ij['p1_int'] and r_j['p2_int'] == r_ij['p2_int']:
                            if r_j['c_int'] != r_ij['c_int']:
                                ACFS_ij.cleavage_reactions.append(r_j)
                            else:
                                break
                #print(i,j)
                #print(ACFS_i.genome_to_catalysts(), ACFS_j.genome_to_catalysts())
                #for r in ACFS_ij.ligand_reactions:
                #    print(r)
                #for r in ACFS_ij.cleavage_reactions:
                #    print(r)
                #print("-----------------------------")
                #for r in ACFS_i.ligand_reactions:
                #    print(r)
                #for r in ACFS_i.cleavage_reactions:
                #    print(r)
                #print("-----------------------------")
                #for r in ACFS_j.ligand_reactions:
                #    print(r)
                #for r in ACFS_j.cleavage_reactions:
                #    print(r)
                # here we solve the kinetic model
                # multiple times -> average different final
                # configurations
                if p.fitness_eval == "compute":
                    ACFS_ij.set_chemical_kinetics_solver()
                #print(ACFS_i.fitness, ACFS_j.fitness)
                # payoff matrix
                self.a_ij[i,j] = ACFS_ij.fitness
                self.a_ij[j,i] = ACFS_ij.fitness
        #
        # write payoff matrix
        f = open(out_file, 'w')
        for i in range(self.size):
            for j in range(self.size):
                f.write("%d     " % i + "%d     " % j + "%.17f\n" % self.a_ij[i,j])
        f.close()
    #
    # extract fitness from file
    def extract_fitness_from_file(self, inp_file, ACF_set):
        # initialize payoff
        self.a_ij = np.zeros((self.size, self.size))
        # open file
        f = open(inp_file, 'r')
        lines = f.readlines()
        # read data
        for line in lines:
            line = line.strip().split()
            i = int(line[0])
            j = int(line[1])
            self.a_ij[i,j] = float(line[2])
        # close file
        f.close()
        # set fitness
        for i in range(self.size):
            ACFS = ACF_set[i]
            ACFS.fitness = self.a_ij[i,i]
    #
    #  compute fitness
    def compute_fitness(self, x_t):
        fitness = np.zeros(len(x_t))
        # compute fitness
        # f(i) = sum_j a_ij x_j
        for i in range(self.size):
            for j in range(self.size):
                fitness[i] += self.a_ij[i,j] * x_t[j]
        return fitness