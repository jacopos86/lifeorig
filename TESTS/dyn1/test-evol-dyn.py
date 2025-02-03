import matplotlib.pyplot as plt
import numpy as np
from lifeorig.fitness_distr import fitness_distr
from lifeorig.set_rndm_matrix import random_matrix
from lifeorig.mutation_rate import zero_mutation

# set parameters

T = 100.
dt = 0.05
n = 30
sig = 1.
r_mut = 0.1
sig_mut = 0.1

# set solver

# compute kernel QSP eq.
def compute_kernel(Q, f, x):
    F = np.zeros(n)
    phi_x = 0.
    for i in range(n):
        phi_x += f[i] * x[i]
    # compute F
    for i in range(n):
        for j in range(n):
            F[i] += Q[i,j] * f[j] * x[j]
        F[i] = F[i] - phi_x * x[i]
    return F

def ODE_solver(y0, Q, f, dt, nt):
    # this routine solves
    # dyi/dt = Qij(t) fj yj - phi(y)yi -> y real n dimensional vector
    # dy/dt = F(y) y
    # Q(t) is n x n matrix
    # using RK4 algorithm
    yt = np.zeros((n,int(nt/2)))
    yt[:,0] = y0[:]
    # iterate over t
    for i in range(int(nt/2)-1):
        y = np.zeros(n)
        y[:] = yt[:,i]
        # K1
        F = compute_kernel(Q[:,:,2*i], f[:,2*i], y)
        K1 = dt * F
        y1 = y + K1 / 2.
        # K2
        F = compute_kernel(Q[:,:,2*i+1], f[:,2*i+1], y1)
        K2 = dt * F
        y2 = y + K2 / 2.
        # K3
        F = compute_kernel(Q[:,:,2*i+1], f[:,2*i+1], y2)
        K3 = dt * F
        y3 = y + K3
        # K4
        F = compute_kernel(Q[:,:,2*i+2], f[:,2*i+2], y3)
        K4 = dt * F
        #
        yt[:,i+1] = y[:] + (K1[:] + 2.*K2[:] + 2.*K3[:] + K4[:]) / 6.
    return yt

# time array
def set_time(dt, T):
    # n. time steps
    nt = int(T / dt)
    time = np.linspace(0., T, nt)
    # dense array
    nt = int(T / (dt/2.))
    time_rk4 = np.linspace(0., T, nt)
    return time, time_rk4

# solve quasi species eq
def solve(x0, fitness_func, mutation_obj):
    time, time_rk4 = set_time(dt, T)
    nt = len(time_rk4)
    # set time evolving fitness
    # fitness (n,2*nt)
    fitness_func.set_constant_fitness_over_time(nt)
    f = fitness_func.fitness_oft
    # set mutation matrix
    # Q (n,n,2*nt)
    mutation_obj.set_constant_mutation_over_time(nt)
    Q = mutation_obj.Q_oft
    # call RK4 method
    xt = ODE_solver(x0, Q, f, dt, nt)
    return time, xt

def normalize_matrix(Q):
    for i in range(n):
        s = sum(Q[i,:])
        Q[i,:] = Q[i,:] / s
    return Q

# set quasi species solver

x0 = np.zeros(n)
for i in range(n):
    x0[i] = np.exp(-(i-n/2)**2 / (2.*sig**2))
x0[:] = x0[:] / np.sum(x0)
print(np.sum(x0))
print(x0)
plt.plot(x0)
plt.show()


# fitness

fitness_func = fitness_distr(n)
fitness_func.fitness[int(n/2)]   = 1.0
fitness_func.fitness[int(n/2)-1] = 0.5
fitness_func.fitness[int(n/2)+1] = 0.5
fitness_func.fitness[int(n/2)-2] = 0.5
fitness_func.fitness[int(n/2)+2] = 0.5
fitness_func.fitness[int(n/2)-3] = 0.25
fitness_func.fitness[int(n/2)+3] = 0.25
print(fitness_func.fitness)
plt.plot(fitness_func.fitness)
plt.show()

# mutation

mutation_obj = zero_mutation(n)
mutation_obj.set_mut_matrix()
#mutation_obj.set_rand_matrix(1, r_mut, sig_mut)
#utation_obj.normalize_matrix()
mutation_obj.Q[int(n/2),int(n/2)+1] = .25
mutation_obj.Q[int(n/2),int(n/2)-1] = .25
mutation_obj.Q[int(n/2),int(n/2)+2] = .2
mutation_obj.Q[int(n/2),int(n/2)-2] = .2
mutation_obj.Q[int(n/2),int(n/2)+3] = .15
mutation_obj.Q[int(n/2),int(n/2)-3] = .15
mutation_obj.Q = normalize_matrix(mutation_obj.Q)

mutation_obj.show_Q_matrix()
for i in range(n):
    print(np.sum(mutation_obj.Q[i,:]))

time, xt = solve(x0, fitness_func, mutation_obj)
print(xt.shape)
plt.plot(xt[:,-1])
plt.show()

