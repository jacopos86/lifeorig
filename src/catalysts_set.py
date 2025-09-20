import matplotlib.pyplot as plt
import numpy as np
from math import exp
from read_input import p
from constants import R
from scipy.stats import truncnorm
import pytest

def reaction_rate(Delta, x0, xgr, sig, E_a, T):
    # E_a in mJ/mol
    N = len(xgr)
    rr = np.zeros(N)
    for i in range(N):
        rr[i] = Delta * exp(-(xgr[i] - x0)**2 / (2*sig**2)) * exp(-E_a / (R*T))
    return rr

#
#  This module builds the catalysts set
#

def build_catalysts_rr_set(size_C, size_C_intern):
    #TODO: run over QSP individuals
    # these are required to compute avg. fitness
    # catalist reaction rate distribution inside ACF set
    E_a = p.ACFS_catalyst_prob_distr["E_a"]
    Delta = p.ACFS_catalyst_prob_distr["Delta"]
    sig = p.ACFS_catalyst_prob_distr["sigma_0"]
    T = p.ACFS_catalyst_prob_distr["temperature"]   # K
    #TODO : use internal catalysts set size 
    xgr_ACFS = np.arange(1, size_C_intern+1, 1)
    x0 = np.random.choice(xgr_ACFS)
    rr_ACFset = reaction_rate(Delta, x0, xgr_ACFS, sig, E_a, T)
    assert(len(rr_ACFset) == size_C_intern)
    rr_ACFset_dict = {"xgr": xgr_ACFS, "rr": rr_ACFset}
    # catalyst reaction rate distribution for ext. catalysts
    xgr_ext = np.arange(size_C_intern+1, size_C+1, 1)
    x0 = np.random.choice(xgr_ext)
    E_a = p.ext_catalyst_prob_distr["E_a"]
    Delta = p.ext_catalyst_prob_distr["Delta"]
    sig = p.ext_catalyst_prob_distr["sigma_0"]
    T = p.ext_catalyst_prob_distr["temperature"]   # K
    #TODO : multiply by external concentration factor
    rr_extset = reaction_rate(Delta, x0, xgr_ext, sig, E_a, T)
    assert(len(rr_extset) == size_C-size_C_intern)
    rr_extset_dict = {"xgr": xgr_ext, "rr": rr_extset}
    # plot distribution
    output_file = p.working_dir + "/rr-distribution.pdf"
    plot_reaction_rate_distr(rr_ACFset_dict, rr_extset_dict, output_file)
    '''
    while len(catalyst_set) < size_C:
        i = random.choice(np.arange(1, size_X+1, 1))
        if i not in catalyst_set:
            catalyst_set.append(i)
    '''
    return rr_ACFset_dict, rr_extset_dict

#
#  set catalysts prob. distr.
#

def set_catalysts_prob_distr(size_C, size_C_intern, ratio_C_ACFset):
    assert(size_C >= size_C_intern)
    p1 = 1./size_C_intern
    if size_C == size_C_intern:
        p2 = 0.
    else:
        p2 = 1./(size_C - size_C_intern)
    prob_distr = np.zeros(size_C)
    for i in range(size_C_intern):
        prob_distr[i] = p1 * ratio_C_ACFset
    for i in range(size_C_intern, size_C):
        prob_distr[i] = p2 * (1. - ratio_C_ACFset)
    print(prob_distr)
    assert(pytest.approx(np.sum(prob_distr), rel=1e-9) == 1.)
    return prob_distr

#
#  build ACFS catalyst distribution
#

def build_catalysts_distr_ACFS(size_X, distr_p):
    mean = distr_p[0]
    std_dev = distr_p[1]
    # Calculate the truncation points for the standard normal distribution
    # (a and b are in terms of standard deviations from the mean)
    min_val = 0
    max_val = size_X
    a = (min_val - mean) / std_dev
    b = (max_val - mean) / std_dev
    # Create a truncated normal distribution object
    truncated_dist = truncnorm(a=a, b=b, loc=mean, scale=std_dev)
    return truncated_dist

#
#  extract list of catalysts from distribution
#

def build_catalysts_list(truncated_dist, N):
    samples = truncated_dist.rvs(N)
    samples = np.round(samples).astype(int)
    return samples

#
#  set up single reaction
#

def set_up_reaction(ratio_C_ACFset, rr_ACFset_dict, rr_extset_dict):
    pass