import numpy as np
# The quasi species solver solves
# the quasi species equation given in input
# 1) a list of fitness values
# 2) a random matrix Q describing the mutations
class QuasiSpeciesSolver():
    def __init__(self, Q, f, n, dt, T):
        # random matrix (n,n,2*nt)
        self.Q = Q
        # fitness (n,2*nt)
        self.f = f
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
    # solve quasi species eq
    def solve(self, x0):
        # call RK4 method
        xt = self.ODE_solver(x0, self.Q, self.f, self.dt)
        # visualize solution
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