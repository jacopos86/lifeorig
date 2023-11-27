import numpy as np
#
# set 2 dimensional world
# LXL sites
# 0 : empty
# 1 : can
#-1 : wall
def display_world(world, L):
    for i in range(L):
        str = ''
        for j in range(L):
            if world[i,j] == 1:
                str += ' * '
            elif world[i,j] == -1:
                str += ' | '
            else:
                str += ' 0 '
        print(str)
L = 15
world = np.zeros((L,L), dtype=int)
nc = 10
x = np.random.randint(1,L-1,nc)
y = np.random.randint(1,L-1,nc)
for i in range(len(x)):
    world[x[i],y[i]] = 1
world[0,:] = -1
world[:,0] = -1
world[L-1,:]= -1
world[:,L-1]= -1
display_world(world, L)