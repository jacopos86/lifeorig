import gillespie
import numpy as np
import matplotlib.pyplot as plt
from lifeorig.logging_module import log
#
#  This module implements
#  the Gillespie algorithm
#  for reaction kinetics
#  Stochastic simulation algorithm :
#  0. Initialize (t=t0, x=x0)
#  1. with system in state (x,t) compute aj(x), a0(x)
#  2. generate tau, j
#  3. next reaction t <- t + tau ; x <- x + vj
#  4. record (x,t) return to 1. 
class chemical_kinetics_solver:
    def __init__(self):
        # initial state list -> population of chemicals
        self.initial_state = []
        # reactions stoichiometry
        self.stoichiometry = []
        # propensities list
        self.propensity = []
    # test function
    def test(self):
        # initial state
        self.initial_state = [290, 10, 0]
        # propensities
        self.propensity = [lambda s, i, r: 2*s*i/300,
                           lambda s, i, r: 0.5*i]
        # stoichiometry
        self.stoichiometry = [[-1, 1, 0],
                              [0, -1, 1]]
        # run simulation
        times, measurements = gillespie.simulate(self.initial_state, self.propensity, self.stoichiometry, duration=15)
        log.debug("\t TEST SUCCESSFUL")
        # plot 
        t = np.array(times)
        state = np.array(measurements)
        plt.plot(t, state[:,0], color='k', linewidth=1.5)
        plt.plot(t, state[:,1], color='b', linewidth=1.5)
        plt.plot(t, state[:,2], color='g', linewidth=1.5)
        plt.plot(t, state[:,0]+state[:,1]+state[:,2], '--', color='k', linewidth=1.5)
        plt.grid()
        plt.show()