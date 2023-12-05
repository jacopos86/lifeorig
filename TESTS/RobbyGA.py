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
L = 20
world = np.zeros((L,L), dtype=int)
nc = 12
x = np.random.randint(1,L-1,nc)
y = np.random.randint(1,L-1,nc)
for i in range(len(x)):
    world[x[i],y[i]] = 1
world[0,:] = -1
world[:,0] = -1
world[L-1,:]= -1
world[:,L-1]= -1
display_world(world, L)
#
# define possible actions
# policy given the actual state
#
# s = (L,S,E,W,N) -> a
# a : move south/pick can/move north/move east/move west/stall
states = []
for L in [-1, 0, 1]:
    for S in [-1,0,1]:
        for E in [-1, 0 ,1]:
            for W in [-1, 0, 1]:
                for N in [-1, 0, 1]:
                    states.append([L,S,E,W,N])
n = len(states)
actions = np.random.randint(6, size=n)
# 0 -> MS
# 1 -> PC
# 2 -> MN
# 3 -> ME
# 4 -> MW
# 5 -> S
