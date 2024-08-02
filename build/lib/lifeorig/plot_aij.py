from lifeorig.parser import parser
import numpy as np
import matplotlib.pyplot as plt

#
#  utility function
#  plot a_ij distribution over grid

args = parser.parse_args()
input_file = args.aij_data[0]
output_file= args.aij_out[0]

#  acquire data from file

f = open(input_file, 'r')
lines = f.readlines()

#  save data

data = []
for line in lines:
    line = line.strip().split()
    
    l = [int(line[0]), int(line[1]), float(line[2])]
    data.append(l)
    
#  a[i,j]

siz = 0

for l in data:
    if l[0] > siz:
        siz = l[0]
siz += 1
a_ij = np.zeros((siz, siz))

#  fill matrix

for l in data:
    i = l[0]
    j = l[1]
    a_ij[i,j] = l[2]
    
#  plot data

plt.matshow(a_ij)
plt.colorbar()
plt.savefig(output_file, format="pdf", bbox_inches="tight")
plt.show()