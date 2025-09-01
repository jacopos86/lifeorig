from parser import parser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

#
#   utility function
#   work as splot_net_fitness
#   but it performs a comparison
#   between different avg. fitness as func. 
#   network size

args = parser.parse_args()
input_files = args.acf_data
output_file = args.acf_data_out[0]
X_size = args.X_size
C_size = args.size_cat

#   plot data

fig, axs = plt.subplots(1, 1, layout='constrained', figsize=(10,10))
color = cm.rainbow(np.linspace(0, 1, len(input_files)))

jf = 0
for fil in input_files:
    
    # open file
    f = open(fil, 'r')
    lines = f.readlines()
    
    net_size = []
    fitness = []
    for line in lines:
        line = line.strip().split()
        net_size.append(int(line[1]))
        fitness.append(float(line[2]))
    
    avg_fit = np.zeros(max(net_size)+1)
    ng_fit = np.zeros(max(net_size)+1)
    for i in range(len(net_size)):
        nets = net_size[i]
        avg_fit[nets] += fitness[i]
        ng_fit[nets] += 1
    for i in range(len(avg_fit)):
        if ng_fit[i] == 0:
            avg_fit[i] = 0
        else:
            avg_fit[i] = avg_fit[i] / ng_fit[i]
    x = np.arange(0, max(net_size)+1)
    
    # make plot

    axs.plot(x, avg_fit, '--', linewidth=2)
    axs.scatter(x, avg_fit, alpha=0.5, label='X size : ' + str(X_size[jf]) + ' - C size : ' + str(C_size[jf]))
    
    f.close()
    jf += 1


axs.set_xlabel('size ACF set')
axs.set_ylabel('network <fitness>')
axs.grid(True)
axs.legend(loc='best')
plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.show()