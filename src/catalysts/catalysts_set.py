import numpy as np
from src.input_data.read_input import p
from src.utilities.logging_module import log
from scipy.stats import truncnorm
import pytest

#
#  This module builds the catalysts set
#

def set_inital_catalyst_set(X_set, rates_params):
    # set catalysts based on distribution

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

def build_catalyst_set(X_set, catalyst_params, normalize=True):
    """
    Select catalyst molecules from X_set based on a length distribution.

    Parameters
    ----------
    X_set : list[Molecule]
        set of all molecules
    catalyst_params : dict
        {
            "distribution": "gaussian",
            "center": float,
            "std": float
        }
    normalize : bool
        whether to normalize probabilities

    Returns
    -------
    Y set: set of all possible catalysts
    """
    dist_type = catalyst_params.get("distribution", "gaussian")
    lengths = np.array([len(mol) for mol in X_set])
    # --- compute weights ---
    if dist_type == "gaussian":
        center = catalyst_params.get("center")
        L0 = lengths[center]
        std = catalyst_params.get("std")
        weights = np.exp(- (lengths - L0)**2 / (2 * std**2))
    else:
        log.error(f"Unknown distribution: {dist_type}")
    # --- normalize ---
    if normalize:
        weights = weights / weights.sum()
    # --- sample exactly n catalysts ---
    indices = np.random.choice(
        len(X_set),
        size=catalyst_params.get("set_size"),
        replace=False,
        p=weights
    )
    catalyst_set = [X_set[i] for i in indices]
    log.info("\n")
    log.info("\t CATALYST SET: \n")
    for X in X_set:
        if X in catalyst_set:
            log.info(f"\t {X.show_sequence()}")
    log.info("\t " + p.sep)
    return catalyst_set

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