import math
from math import sin
from random import random, seed
def f(x):
    return (sin(1/(x*(2-x))))**2
# seed(123)
N = 10000
count = 0
for i in range(N):
    x = 2*random()
    y = random()
    if y<f(x):
        count += 1
I = 2*count/N
# print(I)

import numpy as np


def total_energy(state, J):  # not periodic
    N = len(state)  # size of lattice
    energy = 0  # initialize energy
    # try drawing this first so that you'll know the values of si_sj (be careful with your indexing)
    for i in np.arange(N):
        for j in np.arange(N):
            if i <= j:
                # case1: 4 corner cells - have 2 adjacent cells
                if (i == 0) & (j == 0):
                    si_sj = state[i, j] * (state[i, j + 1] + state[i + 1, j])
                elif (i == 0) & (j == N - 1):
                    si_sj = state[i, j] * (state[i, j - 1] + state[i + 1, j])
                elif (i == N - 1) & (j == 0):
                    si_sj = state[i, j] * (state[i, j + 1] + state[i - 1, j])
                elif (i == N - 1) & (j == N - 1):
                    si_sj = state[i, j] * (state[i, j - 1] + state[i - 1, j])
                # case2: side cells - have 3 adjacent cells
                elif i == 0:
                    si_sj = state[i, j] * (state[i, j + 1] + state[i, j - 1] + state[i + 1, j])
                elif j == N - 1:
                    si_sj = state[i, j] * (state[i, j - 1] + state[i + 1, j] + state[i - 1, j])
                elif i == N - 1:
                    si_sj = state[i, j] * (state[i, j + 1] + state[i, j - 1] + state[i - 1, j])
                elif j == 0:
                    si_sj = state[i, j] * (state[i, j + 1] + state[i + 1, j] + state[i - 1, j])
                # case3: inner cells - have 4 adjacent cells
                else:
                    si_sj = state[i, j] * (state[i, j + 1] + state[i + 1, j] + state[i - 1, j] + state[i, j - 1])
                energy += -1. * si_sj
    return J * energy


def delta_E(y, x, state):
    E_i = total_energy(state, J)  # calculate energy before flipping

    state[y, x] = -1 * state[y, x]  # flip the chosen random state
    E_j = total_energy(state, J)  # calculate new energy after flipping

    return E_j - E_i  # return difference of the two energies


from numpy.random import random
from math import exp


def metropolis_flip(change_E, y, x, state):
    global k_B
    global T
    beta = 1/(k_B*T)  # define beta
    if random() < exp(-beta * change_E):
        state[y, x] = -1*state[y, x]  # accept the flip
    return state


import numpy as np


def magnetization(state):  # calculate the total magnetization of the system
    # global N
    return np.sum(state)


from numpy import sum, zeros
from numpy.random import choice, randint

N = 20  # size of lattice

J = 1  # interaction constant
k_B = 1  # Boltzmann constant
T = 1  # Temperature

seed(123)
state = choice([-1, 1], size=(N, N))  # define NxN lattice which contains 1 or -1 (hint: use numpy.random package)

print(sum(state))
iterations = 1000  # set number of iterations required (use small values to start)
M = []
for i in range(iterations + 1):
    y, x = randint(N), randint(N)  # choose a random spin from the lattice
    change_E = delta_E(y, x, state)  # calculate the change in energy
    state = metropolis_flip(change_E, y, x, state)  # decide whether to flip the state or not

    total_magnetization = magnetization(state)  # compute magnetization of the whole lattice
    M.append(total_magnetization)
import matplotlib.pyplot as plt

# this is for 1000 iterations only, you need to do 1,000,000 ata
# plot magnetization
plt.plot(M)
plt.show()

# help(np.random.choice)
