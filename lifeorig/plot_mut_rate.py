from lifeorig.parser import parser
import matplotlib.pyplot as plt

#
#   utility function
#   plotting the mutation rate

args = parser.parse_args()
input_file = args.mut_inp[0]
output_file= args.mut_out[0]

#   open file

f = open(input_file, 'r')
lines = f.readlines()
t = []
mu_oft = []

for line in lines:
    line = line.strip().split()
    t.append(float(line[0]))
    mu_oft.append(float(line[1]))
    
f.close()

#  plot data

fig, axs = plt.subplots(1, 1, layout='constrained')
axs.plot(t, mu_oft, '--', c='k', linewidth=2)
axs.grid(True)
axs.set_xlabel('t')
axs.set_ylabel('$\mathrm{\mu(t)}$')

#plt.legend(ncol=2, loc='best')
plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.show()