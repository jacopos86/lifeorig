from lifeorig.parser import parser
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm

#
#  utility function
#  plot ACF_data.txt distribution
#  at a given ACFS size

args = parser.parse_args()
input_files = args.acf_data
output_file = args.acf_data_out
acf_size = args.size_acfs
for i in range(len(acf_size)):
    acf_size[i] = int(acf_size[i])
cat_size = args.size_cat

#  acquire data from file

avg_fit = []
fit_acfs = np.zeros((len(input_files), len(acf_size)))

jf = 0
for fil in input_files:
    
    # open file
    f = open(fil, 'r')
    lines = f.readlines()
    
    x = 0.
    n = 0
    n_acfs = np.zeros(len(acf_size), dtype=int)
    f_acfs = np.zeros(len(acf_size))
    for line in lines:
        line = line.strip().split()
        x += float(line[2])
        siz= int(line[1])
        n += 1
        if siz in acf_size:
            ii = acf_size.index(siz)
            n_acfs[ii] += 1
            f_acfs[ii] += float(line[2])
    avg_fit.append(x/n)
    for j in range(len(acf_size)):
        if n_acfs[j] > 0:
            fit_acfs[jf,j] = f_acfs[j] / n_acfs[j]
    jf += 1
    f.close()
    
# plot data

out_file = output_file[0]

fig, axs = plt.subplots(1, 1, layout='constrained')
axs.scatter(cat_size, avg_fit, alpha=0.5)
axs.plot(cat_size, avg_fit, '--', c='k', linewidth=2)
color = cm.rainbow(np.linspace(0, 1, len(acf_size)))
for j in range(len(acf_size)):
    axs.scatter(cat_size, fit_acfs[:,j], alpha=0.1)
    axs.plot(cat_size, fit_acfs[:,j], '-', c=color[j], linewidth=1, label='ACFS size : '+str(acf_size[j]))
axs.grid(True)
axs.set_xlabel('size catalysts set')
axs.set_ylabel('<fitness>')

plt.legend(ncol=2, loc='best')
plt.savefig(out_file, format="pdf", bbox_inches="tight")
plt.show()

#  plot the fitness distribution

out_file = output_file[1]
fig, axs = plt.subplots(len(input_files), 1, layout='constrained', figsize=(10,10))
color = cm.rainbow(np.linspace(0, 1, len(input_files)))

jf = 0
for fil in input_files:
    
    # open file
    f = open(fil, 'r')
    lines = f.readlines()
    
    fitn = []
    for line in lines:
        line = line.strip().split()
        fitn.append(float(line[2]))
    
    # make hist
    
    axs[jf].hist(fitn, bins=1000, color=color[jf], label='catalysts set size : '+str(cat_size[jf]))
    
    # set final plot
    
    axs[jf].grid(True)
    axs[jf].set_xlabel('fitness')
    axs[jf].set_xlim([0.025, 1.])
    axs[jf].set_ylim([0.,90])
    axs[jf].legend(loc='best')
    
    f.close()
    jf += 1

plt.savefig(out_file, format="pdf", bbox_inches="tight")
plt.show()