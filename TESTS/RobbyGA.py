import numpy as np
import sys
from random import randint
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
nc = 100
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
# s = (C,S,E,W,N) -> a
# a : move south(MS)/move north(MN)/move east(ME)/move west(MW)/
#     stay put(SP)/try pick can (TPC)/move random (MR)
states = []
for C in [-1, 0, 1]:
    for S in [-1,0,1]:
        for E in [-1, 0 ,1]:
            for W in [-1, 0, 1]:
                for N in [-1, 0, 1]:
                    states.append([C,S,E,W,N])
n = len(states)
print("possible states: " + str(n))
actions = np.random.randint(7, size=n)
# 0 -> MS
# 1 -> MN
# 2 -> ME
# 3 -> MW
# 4 -> SP
# 5 -> TPC
# 6 -> MR
# strategy : map between state space and actions
def strategy_fitness(x0, maxiter):
    # compute strategy fitness
    # x0 : initial robot location
    x = x0
    iter = 0
    reward = 0
    while (iter < maxiter):
        state = []
        C = world[x[0],x[1]]
        state.append(C)
        S = world[x[0],x[1]-1]
        state.append(S)
        E = world[x[0]+1,x[1]]
        state.append(E)
        W = world[x[0]-1,x[1]]
        state.append(W)
        N = world[x[0],x[1]+1]
        state.append(N)
        # state index
        ist = (N+1) + (W+1)*3 + (E+1)*9 + (S+1)*27 + (C+1)*81
        r, y = update_world_location(ist, x)
        reward += r
        if r is None:
            print("action not recognized")
            sys.exit(1)
        x = y
        iter += 1
    return reward
# update robot state
def update_world_location(ist, x):
    action = actions[ist]
    if action == 6:
        # MR -> move randomly
        action = randint(0, 4)
    r = None
    y = None
    #print("action= ", action)
    # perform action
    if action == 0:
        # MS -> move south
        if x[1] == 1:
            r = -5
            y = x
        else:
            r = 0
            y = [x[0], x[1]-1]
    if action == 1:
        # MN -> move north
        if x[1] == L-2:
            r = -5
            y = x
        else:
            r = 0
            y = [x[0], x[1]+1]
    if action == 2:
        # ME -> move east
        if x[0] == L-2:
            r = -5
            y = x
        else:
            r = 0
            y = [x[0]+1, x[1]]
    if action == 3:
        # MW -> move west
        if x[0] == 1:
            r = -5
            y = x
        else:
            r = 0
            y = [x[0]-1, x[1]]
    if action == 4:
        # SP -> stay put
        r = 0
        y = x
    if action == 5:
        # TPC -> try pick can
        if world[x[0],x[1]] == 0:
            r = -1
            y = x
        elif world[x[0],x[1]] == 1:
            r = 10
            y = x
            world[x[0],x[1]] = 0
        else:
            print("pick can on the wall")
            sys.exit(1)
    return r, y
#
# compute strategy fitness
x0 = [randint(1,L-2), randint(1,L-2)]
fitness = strategy_fitness(x0, 100)
print(fitness)
display_world(world, L)