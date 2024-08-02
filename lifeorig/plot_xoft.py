from lifeorig.parser import parser
import matplotlib.pyplot as plt
from lifeorig.logging_module import log
import numpy as np
from matplotlib.pyplot import cm

#
#  utility function
#  plot the ACF_data.txt
#  size food set cycle vs fitness

args = parser.parse_args()
input_file = args.xt_data[0]
output_file= args.xt_data_out[0]

# arrays

x_oft = []

# open file

f = open(input_file, 'r')
data = f.readlines()
line = data[0].strip().split()
n = len(line)
T = len(data)

# x(t)

x_oft= np.zeros((n, T))
t = 0
for line in data:
    line = line.strip().split()
    j = 0
    for x in line:
        x_oft[j,t] = float(x)
        j += 1
    t += 1
    
# plot data

fig, axs = plt.subplots(1, 1, layout='constrained', figsize=(6,6))
color = cm.rainbow(np.linspace(0, 1, n))
t = np.arange(0, T, 1)

# plot x(t)

log.info("\t T simulations: " + str(len(data)))

for j in range(n-1):
    
    # make plot

    axs.plot(t, x_oft[j,:], '-', c=color[j], linewidth=2, alpha=0.8)

# make plot

axs.plot(t, x_oft[-1,:], '--', c='k', linewidth=2, alpha=0.8)

# plot data

axs.set_xlabel('T')
axs.set_ylabel('x(T)')
axs.grid(True)

plt.savefig(output_file, format="jpeg", bbox_inches="tight")
plt.show()