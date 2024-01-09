import numpy as np
import matplotlib.pyplot as plt
# set parameters
sigma = 1.        # std. dev.
mu = 10.          # mean
tau = 0.05        # time const.
dt = 0.001        # time step
T = 1.            # total time
n = int(T/dt)     # n. time steps
t = np.linspace(0., T, n)
#
sigma_bis = sigma * np.sqrt(2./tau)
sqrt_dt = np.sqrt(dt)
x = np.zeros(n)
#
for i in range(n-1):
    x[i+1] = x[i] + dt * (-(x[i] - mu) / tau) + \
        sigma_bis * sqrt_dt * np.random.randn()
fig, ax = plt.subplots(1, 1, figsize=(8,4))
ax.plot(t, x, lw=2)
plt.show()