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
    def build_X_set(self, size_X):
        self.X_set = list(np.arange(1, size_X+1, 1))
        log.info("\t X set : " + str(self.X_set))
    def set_stoichiometry(self, reaction_set, size_X):
        for r in reaction_set:
            stch = np.zeros(size_X, dtype=int)
            if 'r1_int' in r and 'r2_int' in r and 'p_int' in r:
                i1 = self.X_set.index(r['r1_int'])
                i2 = self.X_set.index(r['r2_int'])
                p  = self.X_set.index(r['p_int'])
                stch[i1]  =-1
                stch[i2] +=-1
                stch[p]  += 1
            elif 'r_int' in r and 'p1_int' in r and 'p2_int' in r:
                i = self.X_set.index(r['r_int'])
                p1= self.X_set.index(r['p1_int'])
                p2= self.X_set.index(r['p2_int'])
                stch[i]  =-1
                stch[p1]+= 1
                stch[p2]+= 1
            self.stoichiometry.append(stch)
        print(len(self.stoichiometry), len(reaction_set))
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