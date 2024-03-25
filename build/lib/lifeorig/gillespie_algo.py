import gillespie
import numpy as np
import matplotlib.pyplot as plt
from lifeorig.logging_module import log
from lifeorig.read_input import p
from matplotlib.pyplot import cm
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
        self.X_mass = []
        for x in self.X_set:
            mx = bin(x).count('1')
            self.X_mass.append(mx)
        log.info("\t X set masses : " + str(self.X_mass)) 
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
    def set_initial_population(self, F_set):
        size_F = len(F_set)
        nF = np.zeros(size_F, dtype=int)
        nF[:] = int(p.n0 / size_F)
        iF = 0
        while (sum(nF) < p.n0):
            nF[iF] += 1
            if iF == size_F-1:
                iF = 0
            else:
                iF+= 1
        log.info("\t initial Food set population : " + str(nF))
        self.initial_state = [None]*len(self.X_set)
        for i in self.X_set:
            ind = self.X_set.index(i)
            if i in F_set:
                self.initial_state[ind] = nF[ind]
            else:
                self.initial_state[ind] = 0
    def set_propensity(self, reaction_set, size_X):
        # run over reactions
        for r in reaction_set:
            if 'r1_int' in r and 'r2_int' in r and 'p_int' in r:
                func = 'lambda '
                r1i = self.X_set.index(r['r1_int'])
                r2i = self.X_set.index(r['r2_int'])
                ci  = self.X_set.index(r['c_int'])
                # prefactor
                for i in range(1, size_X):
                    func += 'x'+str(i)+', '
                func += 'x'+str(size_X)+': '
                if r1i == r2i:
                    func += '(x'+str(r1i+1)+'-1)*x'+str(r1i+1)+'*x'+str(ci+1)+'/2'
                else:
                    func += 'x'+str(r1i+1)+'*x'+str(r2i+1)+'*x'+str(ci+1)+'/2'
                self.propensity.append(eval(func))
            if 'r_int' in r and 'p1_int' in r and 'p2_int' in r:
                func = 'lambda x: '
                ri = self.X_set.index(r['r_int'])
                ci = self.X_set.index(r['c_int'])
                func = 'lambda '
                for i in range(1, size_X):
                    func += 'x'+str(i)+', '
                func += 'x'+str(size_X)+': x'+str(ri+1)+'*x'+str(ci+1)
                self.propensity.append(eval(func))
            #input = tuple(self.initial_state)
            #p = self.propensity[-1]
    def solve(self):
        # run simulation
        log.info("\t " + p.sep)
        log.info("\t START KINETIC MODEL SIMULATION")
        state_t = [None]*p.nconfig
        # run over config.
        for ic in range(p.nconfig):
            times, measurements = gillespie.simulate(self.initial_state, self.propensity, self.stoichiometry, duration=25)
            if ic == 0:
                self.t = np.array(times)
            state_t[ic] = np.array(measurements)
        self.avg_state_t = np.zeros(state_t[0].shape)
        for ic in range(p.nconfig):
            self.avg_state_t += state_t[ic]
        self.avg_state_t = self.avg_state_t / p.nconfig
        log.info("\t END KINETIC SIMULATION")
        log.info("\t " + p.sep)
    def show(self, target_molecules):
        # plot state
        n = len(target_molecules)
        color = iter(cm.rainbow(np.linspace(0, 1, n)))
        for i in range(n):
            ml = target_molecules[i]
            c = next(color)
            plt.plot(self.t, self.X_mass[ml]*self.avg_state_t[:,ml], color=c, linewidth=1.5, label=str(ml))
        # n. molecules
        sm_t = np.zeros(len(self.t))
        for i in range(len(self.t)):
            sm_t[i] = np.sum(self.avg_state_t[i,:])
        plt.plot(self.t, sm_t, '--', color='k', linewidth=2, label='n. molecules')
        # molecules total mass
        size_X = len(self.X_set)
        mm_t = np.zeros(len(self.t))
        for i in range(len(self.t)):
            for j in range(size_X):
                mm_t[i] += self.X_mass[j] * self.avg_state_t[i,j]
        plt.plot(self.t, mm_t, '--', color='b', linewidth=2, label='total mass')
        plt.grid()
        plt.legend()
        plt.show()
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