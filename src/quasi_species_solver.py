import numpy as np
from src.read_input import p
from src.mutation_rate import compute_hamm_dist_matrix
#
# The quasi species solver solves
# the quasi species equation given in input
# 1) a list of fitness values
# 2) a random matrix Q describing the mutations
#
#  function to return
#  the required type of quasi species
#  solver
def BuildQuasiSpeciesSolver(n, dt, T):
    if p.EvolutionaryGameDyn:
        return QuasiSpeciesSolverEvolGDyn(n, dt, T)
    else:
        return QuasiSpeciesSolver(n, dt, T)
    
#
#  base class
class QuasiSpeciesBaseClass(object):
    def __init__(self, n, dt, T):
        # n. of chemistries
        self.n = n
        # time
        self.dt = dt
        self.set_time(dt, T)
    # time array
    def set_time(self, dt, T):
        # n. time steps
        nt = int(T / dt)
        self.time = np.linspace(0., T, nt)
        # dense array
        nt = int(T / (dt/2.))
        self.time_rk4 = np.linspace(0., T, nt)
    #
    # save data on file
    def save_data(self, xt, file_name, mut_file_name):
        # file name
        f = open(file_name, 'w')
        # write data file
        nt = len(self.time_rk4)
        for t in range(int(nt/2)):
            for j in range(self.n):
                f.write("%.17f               " % xt[j,t])
            f.write("%.17f" % sum(xt[:,t]))
            f.write("\n")
        f.close()
        # mutation matrix
        f = open(mut_file_name, 'w')
        # write data file
        nt = len(self.time_rk4)
        for t in range(int(nt/2)):
            for i in range(self.n):
                for j in range(i, self.n):
                    f.write("%.17f               " % self.Q[i,j,t])
            f.write("\n")
        f.close()

#
#  standard Quasi species solver
#  class
#
class QuasiSpeciesSolver(QuasiSpeciesBaseClass):
    def __init__(self, n, dt, T):
        super(QuasiSpeciesSolver, self).__init__(n, dt, T)
    # solve quasi species eq
    def solve(self, x0, fitness_func, mutation_obj, ACF_distr):
        nt = len(self.time_rk4)
        # set time evolving fitness
        # fitness (n,2*nt)
        fitness_func.set_constant_fitness_over_time(nt)
        self.f = fitness_func.fitness_oft
        # set mutation matrix
        # Q (n,n,2*nt)
        if p.mutation_typ == "zero":
            mutation_obj.set_constant_mutation_over_time(nt)
        elif p.mutation_typ == "dist":
            n = len(ACF_distr[0].genome)
            # compute Hamming distance
            HDij = compute_hamm_dist_matrix(ACF_distr)
            n = min(n, HDij.max())
            mu_oft = mutation_obj.set_site_mut_oft(p.r_mut, p.A_mut, p.w_mut, p.tau_mut, self.time_rk4)
            mutation_obj.set_mut_matrix_oft(mu_oft, n, HDij, nt)
            out_file = p.working_dir + "/mu_oft.txt"
            mutation_obj.write_muoft_to_file(self.time_rk4, mu_oft, out_file)
        elif p.mutation_typ == "random":
            mutation_obj.set_constant_mutation_over_time(nt)
        self.Q = mutation_obj.Q_oft
        # call RK4 method
        xt = self.ODE_solver(x0, self.Q, self.f, self.dt)
        return xt
    # compute kernel QSP eq.
    def compute_kernel(self, Q, f, x):
        F = np.zeros(self.n)
        phi_x = 0.
        for i in range(self.n):
            phi_x += f[i] * x[i]
        # compute F
        for i in range(self.n):
            for j in range(self.n):
                F[i] += Q[i,j] * f[j] * x[j]
            F[i] = F[i] - phi_x * x[i]
        return F
    # solver for the 
    # quasi species diff. eq.
    def ODE_solver(self, y0, Q, f, dt):
	    # this routine solves
	    # dyi/dt = Qij(t) fj yj - phi(y)yi -> y real n dimensional vector
        # dy/dt = F(y) y
	    # Q(t) is n x n matrix
        # using RK4 algorithm
        nt = len(self.time_rk4)
        yt = np.zeros((self.n,int(nt/2)))
        yt[:,0] = y0[:]
        # iterate over t
        for i in range(int(nt/2)-1):
            y = np.zeros(self.n)
            y[:] = yt[:,i]
            # K1
            F = self.compute_kernel(Q[:,:,2*i], f[:,2*i], y)
            K1 = dt * F
            y1 = y + K1 / 2.
            # K2
            F = self.compute_kernel(Q[:,:,2*i+1], f[:,2*i+1], y1)
            K2 = dt * F
            y2 = y + K2 / 2.
            # K3
            F = self.compute_kernel(Q[:,:,2*i+1], f[:,2*i+1], y2)
            K3 = dt * F
            y3 = y + K3
            # K4
            F = self.compute_kernel(Q[:,:,2*i+2], f[:,2*i+2], y3)
            K4 = dt * F
            #
            yt[:,i+1] = y[:] + (K1[:] + 2.*K2[:] + 2.*K3[:] + K4[:]) / 6.
        return yt
    
#
#  Quasi species solver class
#  Evolutionary Game dyn.
#
class QuasiSpeciesSolverEvolGDyn(QuasiSpeciesBaseClass):
    def __init__(self, n, dt, T):
        super(QuasiSpeciesSolverEvolGDyn, self).__init__(n, dt, T)
    # solve Q species eq.
    def solve(self, x0, fitness_func, mutation_obj, ACF_distr):
        nt = len(self.time_rk4)
        # time dependent fitness
        # must be computed run time
        # set the mutation matrix
        # Q (n,n,2*nt)
        if p.mutation_typ == "zero":
            mutation_obj.set_constant_mutation_over_time(nt)
        elif p.mutation_typ == "dist":
            n = len(ACF_distr[0].genome)
            # compute Hamming distance
            HDij = compute_hamm_dist_matrix(ACF_distr)
            n = min(n, HDij.max())
            mu_oft = mutation_obj.set_site_mut_oft(p.r_mut, p.A_mut, p.w_mut, p.tau_mut, self.time_rk4)
            mutation_obj.set_mut_matrix_oft(mu_oft, n, HDij, nt)
            out_file = p.working_dir + "/mu_oft.txt"
            mutation_obj.write_muoft_to_file(self.time_rk4, mu_oft, out_file)
        elif p.mutation_typ == "random":
            mutation_obj.set_constant_mutation_over_time(nt)
        self.Q = mutation_obj.Q_oft
        # call RK4 method
        xt = self.ODE_solver(x0, self.Q, fitness_func, self.dt)
        return xt
    # compute kernel QSP eq.
    def compute_kernel(self, Q, f, x):
        F = np.zeros(self.n)
        phi_x = 0.
        for i in range(self.n):
            phi_x += f[i] * x[i]
        # compute F
        for i in range(self.n):
            for j in range(self.n):
                F[i] += Q[i,j] * f[j] * x[j]
            F[i] = F[i] - phi_x * x[i]
        return F
    # solver for the 
    # quasi species diff. eq.
    def ODE_solver(self, y0, Q, fit_func, dt):
	    # this routine solves
	    # dyi/dt = Qij(t) fj yj - phi(y)yi -> y real n dimensional vector
        # dy/dt = F(y) y
	    # Q(t) is n x n matrix
        # using RK4 algorithm
        nt = len(self.time_rk4)
        yt = np.zeros((self.n,int(nt/2)))
        yt[:,0] = y0[:]
        # iterate over t
        for i in range(int(nt/2)-1):
            y = np.zeros(self.n)
            y[:] = yt[:,i]
            # K1
            f = fit_func.compute_fitness(y)
            F = self.compute_kernel(Q[:,:,2*i], f, y)
            K1 = dt * F
            y1 = y + K1 / 2.
            # K2
            f = fit_func.compute_fitness(y1)
            F = self.compute_kernel(Q[:,:,2*i+1], f, y1)
            K2 = dt * F
            y2 = y + K2 / 2.
            # K3
            f = fit_func.compute_fitness(y2)
            F = self.compute_kernel(Q[:,:,2*i+1], f, y2)
            K3 = dt * F
            y3 = y + K3
            # K4
            f = fit_func.compute_fitness(y3)
            F = self.compute_kernel(Q[:,:,2*i+2], f, y3)
            K4 = dt * F
            #
            yt[:,i+1] = y[:] + (K1[:] + 2.*K2[:] + 2.*K3[:] + K4[:]) / 6.
        return yt