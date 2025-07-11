import matplotlib.pyplot as plt
import numpy as np
from math import exp
from src.read_input import p
from src.constants import R

def reaction_rate(Delta, x0, xgr, sig, E_a, T):
    # E_a in mJ/mol
    N = len(xgr)
    rr = np.zeros(N)
    print(E_a/(R*T))
    for i in range(N):
        rr[i] = Delta * exp(-(xgr[i] - x0)**2 / (2*sig**2)) * exp(-E_a / (R*T))
    plt.scatter(xgr, rr)
    plt.show()


#
#  This module builds the catalysts set
#

def build_catalysts_set(size_X, size_C, size_C_ACFset, ratio_C_ACFset):
    catalyst_set = []
    E_a = p.catalyst_prob_distr["E_a"]
    Delta = p.catalyst_prob_distr["Delta"]
    sig = p.catalyst_prob_distr["sigma_0"]
    xgr = np.arange(1, size_C_ACFset+1)
    x0 = 4
    T = p.catalyst_prob_distr["temperature"]   # K
    reaction_rate(Delta, x0, xgr, sig, E_a, T)
    exit()
    # use gaussian prob. distr.
    mean = 10  # Mean of the distribution
    sig  = 2  # Standard deviation of the distribution
    num_samples = 100  # Number of samples to generate
    continuous_samples = np.random.normal(loc=mean, scale=std_dev, size=num_samples)
    while len(catalyst_set) < size_C:
        i = random.choice(np.arange(1, size_X+1, 1))
        if i not in catalyst_set:
            catalyst_set.append(i)
    return catalyst_set