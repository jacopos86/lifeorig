from logging_module import log
from random import randrange
import logging
import os
from mutation_rate import compute_hamming_distance
from read_input import p
from reaction_network import reaction_net_class
from catalysts_set import build_catalysts_distr_ACFS
from plot_catal_distr import plot_ACFS_distr
# This module builds the sample space
# for the evolutionary dynamics
# given a list of network types - it builds the correspondent sample
# size of sample space given in input
#
def check_input_data(size, size_bpol, size_C):
    ntyp = len(size_bpol)
    check = (len(size_C) == ntyp)
    n3 = 0
    for it in range(ntyp):
        n1 = size_bpol[it][1]
        n2 = len(size_C[it])
        check = check and (n1 == n2)
        n3 += n1
    check = check and (size == n3)
    return check

def build_ACFS_networks(size, size_bpol, size_F, size_C):
    # check input data
    assert(check_input_data(size, size_bpol, size_C) == True)
    # build ACFS
    ntyp = len(size_bpol)
    for ityp in range(ntyp):
        net_num = size_bpol[ityp][1]
        net_index = 0
        while net_index < net_num:
            ACFS = reaction_net_class(size_bpol[ityp][0], size_F, size_C[ityp][net_index])
            # set up catalyst set
            size_X = ACFS.size_X
            log.info("\t size_X: " + str(size_X) +
                      " -- catalyst avg: " + str(size_C[ityp][net_index][0]) + 
                      " -- catalyst sig: " + str(size_C[ityp][net_index][1]))
            catalyst_distr = build_catalysts_distr_ACFS(size_X, size_C[ityp][net_index])
            plot_ACFS_distr(catalyst_distr)
            net_index += 1
    exit()
    ACFS.set_binary_polymer_model(catalyst_set)
    ACFS.find_ACF_subset()
    # prepare network plot
    if log.level == logging.INFO:
        file_name = p.working_dir+'/'+str(iter)+'/acfs.html'
        # write food net
        isExist = os.path.isfile(file_name)
        if not isExist:
            os.mkdir(p.working_dir+'/'+str(iter))
        ACFS.show_network(file_name)
    # produce network genome
    ACFS.set_network_genome()
    # here we solve the kinetic model
    # multiple times -> average different final
    # configurations
    ACFS.set_chemical_kinetics_solver()
    log.info("\t SIZE SAMPLE SPACE : " + str(size))
    log.info("\t fitness initial network : " + str(ACFS.fitness))
    # divide size into 10 groups
    ngr_sample = []
    for ig in range(10):
        ngr_sample.append(int(size/10))
    ig = 0
    while (sum(ngr_sample) < size):
        ngr_sample[ig] += 1
        ig += 1
    genome = ACFS.genome
    size_gen = len(genome)
    # sample set
    sample_set = [ACFS]
    # file name
    file_name = p.working_dir + '/'+str(iter) + "/ACF_data.txt"
    f = open(file_name, 'a')
    # append to file
    f.write("%s          " % ACFS.genome + "%d          " % len(ACFS.ACF_set) + "%.7f\n" % ACFS.fitness)
    # run over mutations
    ig = 0
    while (ig < 10):
        for iig in range(ngr_sample[ig]):
            imut = []
            i = 0
            for i in range(ig+1):
                imut.append(randrange(size_gen))
                i += 1
            new_genome = ""
            for ic in range(size_gen):
                int_c = int(genome[ic])
                if ic in imut:
                    new_genome += str((int_c+1)%2)
                else:
                    new_genome += genome[ic]
            # check genome zeros
            ACFS_new = reaction_net_class(size_bpol, size_F, size_C)
            ACFS_new.set_genome_input(new_genome)
            log.info("\t Hamming distance : " + str(compute_hamming_distance(genome, ACFS_new.genome)))
            ACFS_new.set_binary_polymer_from_genome(catalyst_set)
            ACFS_new.find_ACF_subset()
            log.info("\t " + genome)
            log.info("\t " + ACFS_new.genome)
            if log.level <= logging.DEBUG:
                for i in range(len(ACFS.ligand_reactions)):
                    r1 = ACFS.ligand_reactions[i]
                    r2 = ACFS_new.ligand_reactions[i]
                    log.debug("\t catalysts ACF / ACF mutation : " + str(r1['c_int']) + " - " + str(r2['c_int']))
                for i in range(len(ACFS.cleavage_reactions)):
                    r1 = ACFS.cleavage_reactions[i]
                    r2 = ACFS_new.cleavage_reactions[i]
                    log.debug("\t catalysts ACF / ACF mutation : " + str(r1['c_int']) + " - " + str(r2['c_int']))
            # here we solve the kinetic model
            # multiple times -> average different final
            # configurations
            ACFS_new.set_chemical_kinetics_solver()
            # append to file
            f.write("%s          " % ACFS_new.genome + "%d          " % len(ACFS_new.ACF_set) + "%.7f\n" % ACFS_new.fitness)
            # append data to file
            # append new ACF set
            sample_set.append(ACFS_new)
        ig += 1
    # close file
    f.close()
    return sample_set