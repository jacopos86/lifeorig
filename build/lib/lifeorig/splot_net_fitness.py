from lifeorig.parser import parser
import matplotlib.pyplot as plt
from lifeorig.logging_module import log
import numpy as np

#
#  utility function
#  plot the ACF_data.txt
#  size food set cycle vs fitness

args = parser.parse_args()
input_file = args.acf_data[0]
output_file= args.acf_data_out[0]

# arrays
net_size = []
fitness = []

# open file

f = open(input_file, 'r')
data = f.readlines()
for line in data:
    line = line.strip().split()
    net_size.append(int(line[1]))
    fitness.append(float(line[2]))
    
# plot data


plt.show()

# plot average fitness

log.info("\t max net. size : " + str(max(net_size)))
log.info("\t n. simulations: " + str(len(net_size)))

avg_fit = np.zeros(max(net_size)+1)
max_fit = np.zeros(max(net_size)+1)
ng_fit = np.zeros(max(net_size)+1)
for i in range(len(net_size)):
    nets = net_size[i]
    avg_fit[nets] += fitness[i]
    ng_fit[nets] += 1
    if fitness[i] > max_fit[nets]:
        max_fit[nets] = fitness[i]
for i in range(len(avg_fit)):
    if ng_fit[i] == 0:
        avg_fit[i] = 0
    else:
        avg_fit[i] = avg_fit[i] / ng_fit[i]
x = np.arange(0, max(net_size)+1)

# plot data

fig, axs = plt.subplots(2, 1, layout='constrained')
axs[0].scatter(net_size, fitness, alpha=0.5)
axs[0].set_xlabel('size ACF set')
axs[0].set_ylabel('network fitness')
axs[0].grid(True)

axs[0].set_xlabel('size ACF set')
axs[1].set_ylabel('network fitness')
axs[1].plot(x, avg_fit, '--', linewidth=2)
axs[1].scatter(x, avg_fit, alpha=0.5, label='avg. fitness')
axs[1].plot(x, max_fit, ':', c='r', linewidth=2)
axs[1].scatter(x, max_fit, alpha=0.5, color='red', label='max. fitness')
axs[1].grid(True)
plt.legend()

plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.show()